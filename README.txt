TO DO

1) constructor for track_manager5() not working in main{} -- done
2) correct function reducedFrequency5 in test_tracking.cpp; currently returning empty tracks_ -- done
3) if-else added to the function evaluateTracking() giving wrong RMS values -- done
4) reducedFrequency is only copying the reduced frames from 10Hz track to 5Hz track and members label_ and serialization_version_, what about the member track_num_ -- done
5) commented out check lines in main, resolve pointer issue -- done
6) check whether or not timestamps_ of frames have been copied correctly from track_manager to track_manager5 and track_manager2 --
7) check expected behaviour of RMS error for 5Hz and 2Hz colored point cloud trackind --


NOTES

test_tracking.cpp
1) Size of the member estimated_velocity of struct TrackResults is 1 less than the number of frames in the track. This is so because the first time object is seen, velocity is not estimated.
2) The parameter gt_velocities of function getGTVelocities also has a size of 1 less than the number of frames in the track to ensure compatibility with the previous point. 