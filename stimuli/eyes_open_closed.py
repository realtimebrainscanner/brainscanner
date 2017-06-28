import sys
import time
import numpy as np
import argparse
from scipy.io import savemat

from psychopy import visual, core, logging, sound, prefs

	
###########################################################
# Settings
###########################################################

frame_rate = 60.

num_trials = 10
time_closed = 15 # in seconds
time_open = 15 # in seconds
jitter_range = np.arange(-5, 5+1)

num_images = 35

###########################################################
# Parse arguments
###########################################################
parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default = 0)
parser.add_argument('--run', type=int, dest='run', default=1)

args = parser.parse_args()
seed, doRun = args.seed, args.run

filebase = time.strftime('%m-%d-%y_%H-%M-%S')  + '_seed%d' % seed

###########################################################
# Generate pattern for open close paradigm
###########################################################
print('Seed: %d' % seed)
np.random.seed(seed)

image_list = np.random.choice(range(num_images), size=num_trials, replace=False)


#generate jiiter
jitter_open = np.random.choice(jitter_range, size=num_trials, replace=True)
jitter_closed = np.random.choice(jitter_range, size=num_trials, replace=True)

# compute time periods for both conditions (open/close)
open_periods = time_open + jitter_open
close_periods = time_closed + jitter_closed

# translate to number of frames
num_frames_open = open_periods*frame_rate
num_frames_closed = close_periods*frame_rate
num_frames_period = num_frames_open +  num_frames_closed

# convert to integers
num_frames_open = num_frames_open.astype('int')
num_frames_closed = num_frames_closed.astype('int')
num_frames_period = num_frames_period.astype('int')

# build label vector (0 = open, 1 = closed)
labels = np.hstack([np.hstack((np.zeros(open_per), np.ones(close_per))) for open_per, close_per in zip(open_periods, close_periods)])

# save labels in matlab format
paradigm_path = './paradigm_' + filebase
savemat(paradigm_path, {'stim': labels})

if not doRun:
	import sys
	sys.exit()


###########################################################
# Prepare log
#########################################################
logfile_path = './log_' + filebase + '.txt'

globalClock = core.Clock()
logging.LogFile(logfile_path, level=logging.EXP, filemode='w')
logging.setDefaultClock(globalClock)
logging.console = True

def log(msg):
	logging.log(level=logging.EXP, msg=msg)
	print(msg)

###########################################################
# Prepare window and stimuli
###########################################################
win=visual.Window([1920,1080], pos=[0,0], fullscr=True, autoLog=False, monitor="DellLaptop")
fixation = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=0.02,  name='fixation', autoLog=False)

num_images = 10
images = [visual.ImageStim(win, image='images/img%d.jpg' % image_list[idx_image]) for idx_image in range(num_images)] 

open_sound = sound.Sound(value='C', secs=0.5, octave=5, sampleRate=44100, bits=16, autoLog=True, loops=0, name='Open')
closed_sound = sound.Sound(value='C', secs=0.5, octave=5, sampleRate=44100, bits=16, autoLog=True, loops=0, name='Close')
sync_sound = sound.Sound(value='C', secs=0.25, octave=6, sampleRate=44100, bits=16, autoLog=True, loops=0, name='Sync')

###########################################################
# Play synchronization sounds
#########################################################
log('Playing sounds for synchronization')
for idx_sync_sound in range(5):
	for frameN in range(int(frame_rate)):
		if frameN == 0:
			sync_sound.play()

		fixation.draw()
		win.flip()


###########################################################
# Run experiment
#########################################################
log('Experiment started')
for trial_idx in range(num_trials):
	trial_clock = core.Clock()

	log('Beginning trial %d' % trial_idx)

	for frameN in range(num_frames_period[trial_idx]):

		# play sounds
		if frameN == 0:
			open_sound.play()
			pass
		elif frameN == num_frames_open[trial_idx]:
			closed_sound.play()

		# visual
		if 0 <= frameN < num_frames_open[trial_idx]: 
			images[trial_idx].draw()

		if num_frames_open[trial_idx] <= frameN:  # present stim for a different subset
			fixation.draw()

	    # switch buffer
		win.flip()

	log('Ended trial %d' % trial_idx)

	win.clearBuffer()
	win.flip()

	trial_time = trial_clock.getTime()

	print('Time for period: %6.5fms with %4.3f frames per seconds on average' % (trial_time * 1000., num_frames_period[trial_idx]/trial_time) )


###########################################################
# Play sounds to indicate that experiment ended
#########################################################
log('Playing sounds for synchronization')
for idx_sync_sound in range(3):
	for frameN in range(int(frame_rate)):
		if frameN == 0:
			sync_sound.play()

		fixation.draw()
		win.flip()


log('Experimented ended')
win.close()
