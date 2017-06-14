/*
 * test_tracking.cpp
 *
 *  Created on: Jul 30, 2014
 *      Author: davheld
 *
 */

#include <string>
#include <cstdio>
#include <sstream>

#include <boost/math/constants/constants.hpp>
#include <boost/make_shared.hpp>

#include <precision_tracking/track_manager_color.h>
#include <precision_tracking/tracker.h>
#include <precision_tracking/high_res_timer.h>
#include <precision_tracking/sensor_specs.h>

using std::string;

namespace {

const double pi = boost::math::constants::pi<double>();

} // namespace

// Structure for storing estimated velocities for each track.
struct TrackResults {
  int track_num;
  std::vector<Eigen::Vector3f> estimated_velocities;
  std::vector<bool> ignore_frame;
};

// Get the ground-truth velocities for a particular track from the text files in gtFolder.
void getGTVelocities(const string& gt_folder, const int track_num,
                     std::vector<double>* gt_velocities) {
  string filename = gt_folder + "/track%dgt.txt";

  std::ostringstream filename_stream;
  filename_stream << gt_folder << "/track" << track_num << "gt.txt";
  filename = filename_stream.str();

  FILE* fid = fopen(filename.c_str(), "r");

  if (fid == NULL) {
    printf("Cannot open file: %s\n", filename.c_str());
    exit(1);
  }

  double velocity;
  while(fscanf(fid, "%lf\n", &velocity) > 0) {
    gt_velocities->push_back(velocity);
  }

  fclose(fid);
}

// Computen statistics (RMS value) of an array consisting of error in estimated velocity magnitude
void computeErrorStatistics(const std::vector<double>& errors) {
  double sum_sq = 0;

  size_t num_frames = errors.size();

  // Compute the root-mean-square error.
  for (size_t i = 0; i < num_frames; ++i) {
    sum_sq += pow(errors[i], 2);
  }

  const double rms_error = sqrt(sum_sq / errors.size());

  printf("RMS error: %lf m/s\n", rms_error);
}

// Compute error in magnitude of estimated velocity frame-wise for every track and then compute overall RMS error
void evaluateTracking(const std::vector<TrackResults>& velocity_estimates,
                      const string& gt_folder,
                      boost::shared_ptr<std::vector<bool> > filter) {

  std::vector<double> errors;

  int total_framenum = -1;

  // Iterate over all tracks.
  for (size_t i = 0; i < velocity_estimates.size(); ++i) {
    TrackResults track_results = velocity_estimates[i]; // velocity (not speed) estimates for i-th track

    int track_num = track_results.track_num; // the track number of the i-th track

    std::vector<double> gt_velocities;
    // getting reduced frequency 5Hz and 2Hz ground truth velocities
    std::vector<double> gt_velocities5;
    std::vector<double> gt_velocities2;

    getGTVelocities(gt_folder, track_num, &gt_velocities); // get ground truth velocities of the i-th track using its track number and then reading it off of the text file for the i-th track in gtFolder

    for (size_t j = 0; j < gt_velocities.size(); ++j) {
      if (j % 2 != 0) {
        gt_velocities5.push_back(gt_velocities[j]);
      }
    }

    for (size_t j = 0; j < gt_velocities.size(); ++j) {
      if (j % 5 == 4) {
        gt_velocities2.push_back(gt_velocities[j]);
      }
    }

    // Checking which frequency are we working with right now
    // If 10Hz, retain gt_velocities as it is
    // If 5Hz, replace it with gt_velocities5 and if 2Hz, replace it with gt_velocities2
    if (track_results.estimated_velocities.size() == gt_velocities5.size()) {
      gt_velocities = gt_velocities5;
      // printf("Working with 5Hz ground truth velocities\n");
    } else if (track_results.estimated_velocities.size() == gt_velocities2.size()) {
      gt_velocities = gt_velocities2;
      // printf("Working with 2Hz ground truth velocities\n");
    }

    int skipped = 0; // stores number of frames skipped for the i-th track

    // j represents the frame number for the i-th track (qualitatively)
    for (size_t j = 0; j < track_results.estimated_velocities.size(); ++j) {
      total_framenum++;

      // ignore_frame tells whether or not a particular frame has to be ignored, must be taking values as either 0 or 1
      if (track_results.ignore_frame[j]) {
        skipped++;
        continue;
      }

      if (filter && !((*filter)[total_framenum])) {
        continue;
      }

      const Eigen::Vector3f& estimated_velocity =
          track_results.estimated_velocities[j];

      const double estimated_velocity_magnitude = estimated_velocity.norm();
      const double gt_velocity_magnitude = gt_velocities[j-skipped];
      const double error = estimated_velocity_magnitude - gt_velocity_magnitude;

      errors.push_back(error);
    }
  }

  computeErrorStatistics(errors);
}

// Filter to only evaluate on objects within a given distance (in meters).
void getWithinDistance(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const double max_distance, std::vector<bool>& filter) {
  const std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> >& tracks =
      track_manager.tracks_;
  for (size_t i = 0; i < tracks.size(); ++i) {
    // Extract frames.
    const boost::shared_ptr<precision_tracking::track_manager_color::Track>& track = tracks[i];
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > frames =
        track->frames_;
    for (size_t j = 1; j < frames.size(); ++j) {
      boost::shared_ptr<precision_tracking::track_manager_color::Frame> frame = frames[j];

      Eigen::Vector3f centroid = frame->getCentroid();
      const double distance = sqrt(pow(centroid(0), 2) + pow(centroid(1), 2));

      if (distance <= max_distance) {
        filter.push_back(true);
      } else {
        filter.push_back(false);
      }
    }
  }
}

// Ingore frames in the back where half of the car was recorded at
// the beginning of a spin and the other half was recorded at the end of a spin.
// Also ignore frames where the time difference is extremely small, essentially
// because the car moved between the end of one spin to the beginning of the
// next.  For such frames, estimating the velocity is prone to errors
// that should ideally be fixed before the track is passed on to the velocity
// estimator.
void find_bad_frames(const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
                     std::vector<TrackResults>* velocity_estimates) {
  const std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> >& tracks =
      track_manager.tracks_;

  // Iterate over all tracks.
  for (size_t i = 0; i < tracks.size(); ++i) {
    // Extract frames.
    const boost::shared_ptr<precision_tracking::track_manager_color::Track>& track = tracks[i];
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > frames =
        track->frames_;

    // Structure for storing estimated velocities for this track.
    TrackResults& track_estimates = (*velocity_estimates)[i];

    bool skip_next = false;
    double prev_angle = 0;
    double prev_time = 0;

    // Iterate over all frames for this track.
    for (size_t j = 0; j < frames.size(); ++j) {
      boost::shared_ptr<precision_tracking::track_manager_color::Frame> frame = frames[j];

      const Eigen::Vector3f& centroid = frame->getCentroid();
      const double angle = atan2(centroid(1), centroid(0));
      const double angle_diff = fabs(angle - prev_angle);

      const double curr_time = frame->timestamp_;
      const double time_diff = curr_time - prev_time;

      prev_angle = angle;

      if (j > 0) {
        if (angle_diff <= 1 || j == 0) {
          if (!skip_next) {
            if (time_diff >= 0.05) {
              //track_estimates.ignore_frame.push_back(false);
            } else {
              track_estimates.ignore_frame[j-1] = true;
            }
          } else {
            track_estimates.ignore_frame[j-1] = true;
          }
          skip_next = false;
        } else {
          track_estimates.ignore_frame[j-1] = true;
          skip_next = true;
          if (j > 1) {
            track_estimates.ignore_frame[j-2] = true;
          }
        }
      }
    }
  }
}

void track(
           const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
           const precision_tracking::Params& params,
           const bool use_precision_tracker,
           const bool do_parallel,
           std::vector<TrackResults>* velocity_estimates) {
  const std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> >& tracks =
      track_manager.tracks_;

  int total_num_frames = 0;
  for (size_t i = 0; i < tracks.size(); ++i) {
    total_num_frames += tracks[i]->frames_.size();
  }

  const int num_threads = do_parallel ? 8 : 1;

  std::vector<precision_tracking::Tracker> trackers;
  for (int i = 0; i < num_threads; ++i) {
    precision_tracking::Tracker tracker(&params);
    if (use_precision_tracker) {
      tracker.setPrecisionTracker(
          boost::make_shared<precision_tracking::PrecisionTracker>(&params));
    }
    trackers.push_back(tracker);
  }

  velocity_estimates->resize(tracks.size());

  std::ostringstream hrt_title_stream;
  hrt_title_stream << "Total time for tracking " << tracks.size() << " objects";
  precision_tracking::HighResTimer hrt(hrt_title_stream.str(),
                                       do_parallel ? CLOCK_REALTIME :
                                                     CLOCK_PROCESS_CPUTIME_ID);
  hrt.start();

  // Iterate over all tracks.
  #pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < tracks.size(); ++i) {

    precision_tracking::Tracker& tracker = trackers[omp_get_thread_num()];

    // Reset the tracker for this new track.
    tracker.clear();

    // Extract frames.
    const boost::shared_ptr<precision_tracking::track_manager_color::Track>& track = tracks[i];
    const std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > frames =
        track->frames_;

    // Structure for storing estimated velocities for this track.
    TrackResults track_estimates;
    track_estimates.track_num = track->track_num_;

    // Iterate over all frames for this track.
    for (size_t j = 0; j < frames.size(); ++j) {
      const boost::shared_ptr<precision_tracking::track_manager_color::Frame> frame = frames[j];

      // Get the sensor resolution.
      double sensor_horizontal_resolution;
      double sensor_vertical_resolution;
      precision_tracking::getSensorResolution(
            frame->getCentroid(), &sensor_horizontal_resolution,
            &sensor_vertical_resolution);

      // Track object.
      Eigen::Vector3f estimated_velocity;
      tracker.addPoints(frame->cloud_, frame->timestamp_,
                         sensor_horizontal_resolution,
                         sensor_vertical_resolution,
                         &estimated_velocity);

      // The first time we see this object, we don't have a velocity yet.
      // After the first time, save the estimated velocity.
      if (j > 0) {
        track_estimates.estimated_velocities.push_back(estimated_velocity);

        // By default, don't ignore any frames.
        track_estimates.ignore_frame.push_back(false);
      }
    }
    (*velocity_estimates)[i] = track_estimates;
  }

  hrt.stop();
  hrt.print();

  const double ms = hrt.getMilliseconds();
  printf("Mean runtime per frame: %lf ms\n", ms / total_num_frames);
}

void trackAndEvaluate(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder,
    const precision_tracking::Params& params,
    const bool use_precision_tracker,
    const bool track_parallel) {
  // Track all objects and store the estimated velocities.
  std::vector<TrackResults> velocity_estimates;
  track(track_manager, params, use_precision_tracker, track_parallel, &velocity_estimates);

  // Find bad frames that we want to ignore.
  find_bad_frames(track_manager, &velocity_estimates);

  // Evaluate the tracking accuracy for all objects in the sensor field of view.
  boost::shared_ptr<std::vector<bool> > empty_filter;
  evaluateTracking(velocity_estimates, gt_folder, empty_filter);

  // Evaluate the tracking accuracy for nearby objects.
  const double max_distance = 5;
  printf("Evaluating only for objects within %lf m:\n", max_distance);
  boost::shared_ptr<std::vector<bool> > filter(new std::vector<bool>);
  getWithinDistance(track_manager, max_distance, *filter);
  evaluateTracking(velocity_estimates, gt_folder, filter);
}

void testKalman(const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
                const string gt_folder) {
  printf("Tracking objects with the centroid-based Kalman filter baseline. "
         "This method is very fast but not very accurate. Please wait...\n");
  precision_tracking::Params params;
  trackAndEvaluate(track_manager, gt_folder, params, false, false);
}

void testPrecisionTracker2D(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with our precision tracker in 2D (single-threaded). "
         "This method is accurate and fairly fast. Compared to the full 3D version, this method uses much less memory "
         "and is much faster, but is slightly less accurate.  Please wait...\n");
  precision_tracking::Params params;
  trackAndEvaluate(track_manager, gt_folder, params, true, false);
}

void testPrecisionTracker2DParallel(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with our precision tracker in 2D in parallel. "
         "This method is accurate and fairly fast. Compared to the full 3D version, this method uses much less memory "
         "and is much faster, but is slightly less accurate.  Please wait...\n");
  precision_tracking::Params params;
  trackAndEvaluate(track_manager, gt_folder, params, true, true);
}

void testPrecisionTracker3D(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with our precision tracker in 3D (single-threaded). "
         "This method is accurate and fairly fast. Compared to the 2D version, this method uses more memory "
         "and is slower, but is more accurate.  Please wait...\n");
  precision_tracking::Params params;
  params.use3D = true;
  trackAndEvaluate(track_manager, gt_folder, params, true, false);
}

void testPrecisionTrackerColor(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with our precision tracker using color (single-threaded). "
         "This method is a bit more accurate than the version without color but is much slower. Please wait (will be slow)...\n");
  precision_tracking::Params params;
  params.useColor = true;
  trackAndEvaluate(track_manager, gt_folder, params, true, false);
}

// Converting track_manager to reduced 5Hz frequency.
void reducedFrequency5(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    precision_tracking::track_manager_color::TrackManagerColor& track_manager5) {
  printf("\nConverting %zu tracks of colored laser point clouds from 10Hz to 5Hz and storing it in an empty list, currently having %zu tracks\n", track_manager.tracks_.size(), track_manager5.tracks_.size());

  // Extract pointer to tracks for both original 10Hz and currently empty 5Hz list of tracks
  std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> > tracks =
      track_manager.tracks_;
  std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> > tracks5 =
      track_manager5.tracks_;

  // Iterate over all tracks.
  for (size_t i = 0; i < tracks.size(); ++i) {
    boost::shared_ptr<precision_tracking::track_manager_color::Track>& track = tracks[i];
    // Extract frames.
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > frames =
        track->frames_;
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > subsframes;
    for (size_t j = 0; j < frames.size(); ++j) {
      if (j % 2 == 0) continue;
      subsframes.push_back(frames[j]);
    }
    boost::shared_ptr<precision_tracking::track_manager_color::Track> track5(new precision_tracking::track_manager_color::Track(track->label_, subsframes));
    int track5_num = track->track_num_;
    track5->track_num_ = track5_num;
    tracks5.push_back(track5);
  }
  track_manager5.tracks_ = tracks5;
}

// Converting track_manager to reduced 2Hz frequency.
void reducedFrequency2(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    precision_tracking::track_manager_color::TrackManagerColor& track_manager2) {
  printf("\nConverting %zu tracks of colored laser point clouds from 10Hz to 2Hz and storing it in an empty list, currently having %zu tracks\n", track_manager.tracks_.size(), track_manager2.tracks_.size());

  // Extract pointer to tracks for both original 10Hz and currently empty 5Hz list of tracks
  std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> > tracks =
      track_manager.tracks_;
  std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Track> > tracks2 =
      track_manager2.tracks_;

  // Iterate over all tracks.
  for (size_t i = 0; i < tracks.size(); ++i) {
    boost::shared_ptr<precision_tracking::track_manager_color::Track>& track = tracks[i];
    // Extract frames.
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > frames =
        track->frames_;
    std::vector< boost::shared_ptr<precision_tracking::track_manager_color::Frame> > subsframes;
    for (size_t j = 0; j < frames.size(); ++j) {
      if (j % 5 == 0) {
        subsframes.push_back(frames[j]);
      }
    }
    boost::shared_ptr<precision_tracking::track_manager_color::Track> track2(new precision_tracking::track_manager_color::Track(track->label_, subsframes));
    int track2_num = track->track_num_;
    track2->track_num_ = track2_num;
    tracks2.push_back(track2);
  }
  track_manager2.tracks_ = tracks2;
}

// Input to the following functions will be track_manager5/track_manager2 and gt_folder as defined in main()
void testPrecisionTrackerColor5Hz(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with our precision tracker using color (single-threaded) and reduced frequency of 5Hz. "
         "This method should also be less accurate than the version with color at 10Hz but is slightly faster. Please wait (will be slow)...\n");
  precision_tracking::Params params;
  params.useColor = true;
  trackAndEvaluate(track_manager, gt_folder, params, true, false);
}

void testPrecisionTrackerColor2Hz(
    const precision_tracking::track_manager_color::TrackManagerColor& track_manager,
    const string gt_folder) {
  printf("\nTracking objects with DH precision tracker using color (single-threaded) and reduced frequency of 2Hz. "
         "This method should be less accurate than the version with color at 10Hz but much faster. Please wait (will be slow)...\n");
  precision_tracking::Params params;
  params.useColor = true;
  trackAndEvaluate(track_manager, gt_folder, params, true, false);
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    printf("Usage: %s tm_file gt_folder\n", argv[0]);
    return (1);
  }

  string color_tm_file = argv[1];
  string gt_folder = argv[2];

  // Load tracks.
  printf("Loading file: %s\n", color_tm_file.c_str());

  // Use the constructor for class TrackManagerColor which takes as input a filename
  // Function definiton for this particular constructor can be found in track_manager_color.cpp searching TrackManagerColor::TrackManagerColor(const string& filename)
  precision_tracking::track_manager_color::TrackManagerColor track_manager(color_tm_file);
  printf("Found %zu tracks\n", track_manager.tracks_.size());

  // Track objects and evaluate the accuracy.
  printf("Tracking objects - please wait...\n\n");

  // Testing the centroid-based Kalman filter baseline method - should be
  // very fast but not very accurate.
  testKalman(track_manager, gt_folder);

  // Testing our precision tracker - should be very accurate and quite fast.
  testPrecisionTracker2D(track_manager, gt_folder);

  // Testing our precision tracker - should be very accurate and quite fast.
  testPrecisionTracker2DParallel(track_manager, gt_folder);

  // Testing our precision tracker - should be very accurate and quite fast.
  testPrecisionTracker3D(track_manager, gt_folder);

  // Testing our precision tracker with color - should be even more accurate
  // but slow.
  testPrecisionTrackerColor(track_manager, gt_folder);

  // // Use default constructor TrackManagerColor() defined in track_manager_color.cpp, which initializes the member tracks_ to an empty vector of the class Track as defined in TrackManaagerColor.h
  // printf("\nInitializing empty track_manager class for storing and handling reduced 5Hz tracks."
  //        "\nShould have no tracks if successfully initialized.");
  precision_tracking::track_manager_color::TrackManagerColor track_manager5;
  // printf("\nFound %zu tracks compatible with %d serialization version", track_manager5.tracks_.size(), track_manager5.serialization_version_);

  // Convert track_manager and gt_folder to reduced 5Hz frequency.
  printf("\n\nConverting tracks from 10Hz to 5Hz");
  reducedFrequency5(track_manager, track_manager5);

  // printf("\nChecking number of reduced tracks, should be the same as before..."
  //        "\nFound %zu tracks", track_manager5.tracks_.size());
  // printf("\n\nChecking number of frames for randomly selected tracks, should be half of the original");
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[1]->frames_.size(), track_manager5.tracks_[1]->frames_.size(), track_manager.tracks_[1]->track_num_, track_manager5.tracks_[1]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[121]->frames_.size(), track_manager5.tracks_[121]->frames_.size(), track_manager.tracks_[121]->track_num_, track_manager5.tracks_[121]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[256]->frames_.size(), track_manager5.tracks_[256]->frames_.size(), track_manager.tracks_[256]->track_num_, track_manager5.tracks_[256]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[398]->frames_.size(), track_manager5.tracks_[398]->frames_.size(), track_manager.tracks_[398]->track_num_, track_manager5.tracks_[398]->track_num_);

  // Testing DH precision tracker with color - at reduced frequency of 5Hz.
  testPrecisionTrackerColor5Hz(track_manager5, gt_folder);

  // printf("\n\nInitializing empty track_manager class for storing and handling reduced 2Hz tracks."
  //        "\nShould have no tracks if successfully initialized.");
  precision_tracking::track_manager_color::TrackManagerColor track_manager2;
  // printf("\nFound %zu tracks compatible with %d serialization version", track_manager2.tracks_.size(), track_manager2.serialization_version_);

  // Convert track_manager and gt_folder to reduced 5Hz frequency.
  printf("\n\nConverting tracks from 10Hz to 2Hz");
  reducedFrequency2(track_manager, track_manager2);

  // printf("\nChecking number of reduced tracks, should be the same as before..."
  //        "\nFound %zu tracks", track_manager2.tracks_.size());

  // printf("\n\nChecking number of frames for randomly selected tracks, should be one-fifth of the original");
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[1]->frames_.size(), track_manager2.tracks_[1]->frames_.size(), track_manager.tracks_[1]->track_num_, track_manager2.tracks_[1]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[121]->frames_.size(), track_manager2.tracks_[121]->frames_.size(), track_manager.tracks_[121]->track_num_, track_manager2.tracks_[121]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d", track_manager.tracks_[256]->frames_.size(), track_manager2.tracks_[256]->frames_.size(), track_manager.tracks_[256]->track_num_, track_manager2.tracks_[256]->track_num_);
  // printf("\nFound %zu frames and %zu reduced frames in track number %d and %d\n", track_manager.tracks_[398]->frames_.size(), track_manager2.tracks_[398]->frames_.size(), track_manager.tracks_[398]->track_num_, track_manager2.tracks_[398]->track_num_);

  // Testing DH precision tracker with color - at reduced frequency of 2Hz.
  testPrecisionTrackerColor2Hz(track_manager2, gt_folder);

  // std::vector<double> gt_velocities1;
  // std::vector<double> gt_velocities2;
  // std::vector<double> gt_velocities3;
  // std::vector<double> gt_velocities4;
  // getGTVelocities(gt_folder, 1, &gt_velocities1);
  // getGTVelocities(gt_folder, 2792, &gt_velocities2);
  // getGTVelocities(gt_folder, 3978, &gt_velocities3);
  // getGTVelocities(gt_folder, 4702, &gt_velocities4);
  // printf("Checking size of gt_velocities for the same frames as picked before = %zu, %zu, %zu, %zu\n", gt_velocities1.size(), gt_velocities2.size(), gt_velocities3.size(), gt_velocities4.size());

  return 0;
}

