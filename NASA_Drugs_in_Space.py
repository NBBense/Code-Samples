# Drugs in Space: A Scientific Payload for 
# Investigations of RNA Folding, Synthetic 
# Protocell Protein Expression, and 
# Pharmaceutical Nano-Precipitation.

# Import packages
import RPi.GPIO as GPIO
import time
import picamera
from numpy.random import randint

# Function to delay GPIO pin activation/deactivation
def wait(start, duration):
    while time.time()  - start < duration:
        time.sleep(1)

# Initialize starting time value
starting_time = time.time()

# Initialize camera and start recording
camera=picamera.PiCamera()
camera.start_recording('microgravity_experiment%d.h264' % randint(0,50000))

# Assign equipment variables to pins
LED = 3
GSK_SLOW_PUMP = 7
GSK_FAST_PUMP = 11
UMN_PUMPS = 16

# Pin intialization
GPIO.cleanup()
GPIO.setmode(GPIO.BOARD)
GPIO.setup(LED, GPIO.OUT)
GPIO.setup(GSK_SLOW_PUMP, GPIO.OUT)
GPIO.setup(GSK_FAST_PUMP, GPIO.OUT)
GPIO.setup(UMN_PUMPS, GPIO.OUT)

# Assign times for pin deactivation
TURN_ON_GSK_SLOW = 10
TURN_ON_GSK_FAST_AND_UMN = 720
TURN_OFF_UMN = 765
TURN_OFF_GSK_SLOW = 1005
TURN_OFF_CAMERA_LEDS = 2000

# Activate GSK slow pump
wait(starting_time, TURN_ON_GSK_SLOW)
GPIO.output(GSK_SLOW_PUMP, 0)

# Activate GSK fast pump and UMN pumps
wait(staring_time, TURN_ON_GSK_FAST_AND_UMN)
GPIO.output(GSK_FAST_PUMP, 1)
GPIO.output(UMN_PUMPS, 1)

# Turn off UMN pumps
wait(starting_time, TURN_OFF_UMN)
GPIO.output(UMN_PUMPS, 0)

# Turn off GSK slow pump
wait(starting_time, TURN_OFF_GSK_SLOW)
GPIO.output(GSK_SLOW_PUMP, 0)
GPIO.output(GSK_FAST_PUMP, 0)

# Turn off LEDs and camera
wait(starting_time, TURN_OFF_CAMERA_LEDS)
GPIO.output(LED, 0)
camera.stop_recording()

######END OF PROGRAM########
