#!/usr/bin/env python
#
# Heart Rate and Electrocardiagram Simulation
# CSE 6730 Modeling and Simulation Project #2
# Dylan Crocker and Zeaid Hasan
# Georgia Institute of Technology
# May 2016
#

# Imports -------------------------------------------------------------------- #

import sys
import time
import Queue
import numpy as np
import ConfigParser
from uuid import uuid4
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Constants ------------------------------------------------------------------ #

ECG_AMP_MAX = +1.2  # Maximum ECG signal amplitude (mV)
ECG_AMP_MIN = -0.4  # Minimum ECG signal amplitude (mV)
ECG_NOISE_AMP = .2  # ECG Noise signal max amplitude (mV)
START_ACTIVITY = 1  # Start activity event flag
END_ACTIVITY = -1   # End activity event flag

# Classes -------------------------------------------------------------------- #


class Human(object):
    """Store characteristics of the human whose heart will be simulated."""

    def __init__(self, verbose=False):
        """Create a human object with all default parameters."""

        # Attributes to set once during creation of the object.
        self.age = 45
        self.gender = "male"
        self.mass = 80
        self.rhr = 60
        self.vo2max = 11.0

        # State variables that change throughout the simulation
        self.hr = 60
        self.cav = 7
        self.intensity = 1
        self.activity = None

        # Control variables
        self.verbose = verbose

    def description(self):
        """A simple text description of the object."""
        template = "The human subject is {} age {}"
        return template.format(self.gender, self.age)

    def change_activity(self, activity):
        self.activity = activity
        self._adjust_activity_level(activity.met)

    def rest(self):
        self.activity = None
        self._adjust_activity_level(1.0)

    def _adjust_activity_level(self, intensity):
        """Mutate the human object for the new activity level.

        :param intensity: Activity intensity in units of METs.
        """
        self.intensity = intensity
        sv = self._calculate_stroke_volume()
        self.cav = self._calculate_cav(intensity)
        self.hr = 3.5*self.mass*intensity/sv/(self.cav/100)
        if self.verbose:
            print("HR = {} beats/min".format(self.hr))

    def _initial_stroke_volume(self):
        """Calculate the human's resting stroke volume."""
        sv = 3.5*self.mass/self.rhr/(self._calculate_cav(1.0)/100)
        if self.verbose:
            print("Initial SV = {} ml".format(sv))
        return sv  # mL

    def _calculate_stroke_volume(self):
        """Calculate updated stroke volume.

        This uses a linear approximation.

        :return: Stroke volume (ml) scaled for activity level.
        """
        max_increase = 0.65  # 65% increase at max
        sv_init = self._initial_stroke_volume()
        if self.intensity/self.vo2max >= 0.6:
            # After 60% VO2max the SV has plateaued
            sv = sv_init*(1 + max_increase)
        elif self.intensity > 1:
            # Model as linear increase up to plateau
            sv = sv_init*(max_increase*(self.intensity - 1)/(0.6*self.vo2max - 1) + 1)
        else:
            # Keep resting SV
            sv = sv_init
        if self.verbose:
            print("Scaled SV = {} ml".format(sv))
        return sv

    def _calculate_cav(self, intensity):
        """Calculate arteriovenous oxygen content difference (Cav).

        :param intensity: Exercise intensity in units of METs.
        """
        cav = 5.72 + 0.1047*(intensity/self.vo2max*100)
        if self.verbose:
            print("Cav = {} ml/100ml".format(cav))
        return cav


class Activity(object):
    """Represents an activity that will change the average heart rate."""

    def __init__(self, uid):
        """Create an activity object with all default values."""
        self.type = "resting"
        self.met = 1
        self.start = 0
        self.duration = 60
        self.uid = uid

    def description(self):
        """A simple text description of the object."""
        template = "Perform {} activity for {} minutes"
        return template.format(self.type, self.duration)


class Event(object):
    """ """

    def __init__(self, event_id, event_obj):
        self.id = event_id
        self.obj = event_obj


class HeartSimulation(object):
    """Simulate average HR based on activity level"""

    def __init__(self, human, fel, verbose=False, visual=False):
        self.human = human
        self.fel = fel
        self.avg_hr = []
        self.verbose = verbose
        self.visual = visual

    def start_activity_event(self, time_stamp, event):
        """Start activity event handler.

        :param time_stamp: Time stamp of the activity start.
        :param event: Event object containing event information.
        """
        activity = event.obj

        if self.verbose:
            if self.human.activity is not None:
                print("End activity {} at time {}".format(self.human.activity.type, time_stamp))
            print("\nStart activity {} at time {}".format(activity.type, time_stamp))

        # Queue an event that will end the activity
        new_event = Event(END_ACTIVITY, activity)
        self.fel.put((time_stamp + activity.duration, new_event))

        old_hr = self.human.hr
        self.human.change_activity(activity)
        if self.verbose:
            print("delta HR = {}".format(self.human.hr - old_hr))

        # Save the change in HR
        self.avg_hr.append((time_stamp, self.human.hr, activity.met))

    def end_activity_event(self, time_stamp, event):
        """End activity event handler.

        :param time_stamp: Time stamp of the activity start.
        :param event: Event object containing event information.
        """
        activity = event.obj

        # Check to see if the activity is still in progress.
        if activity.uid == self.human.activity.uid:
            if self.verbose:
                print("End activity {} at time {}".format(activity.type, time_stamp))
            old_hr = self.human.hr
            self.human.rest()  # Put the heart back at rest
            if self.verbose:
                print("delta HR = {}".format(self.human.hr - old_hr))

            # Save the change in HR
            self.avg_hr.append((time_stamp, self.human.hr, 1.0))

    def run_simulation(self, output_file_path=None):
        """Run the discrete event heart rate simulation."""

        # Process queued events
        while not self.fel.empty():

            # Get the next event with lowest time stamp value
            now, event = self.fel.get()

            # Call event handlers
            if event.id == START_ACTIVITY:
                self.start_activity_event(now, event)
            elif event.id == END_ACTIVITY:
                self.end_activity_event(now, event)

        # Process the HR data to include transitions.
        # Approximate transitions as 2 min linear transitions discretely stepped every 10 seconds.
        # Assume the HR starts at rest
        temp_events = []
        t_step = 1./6
        prev_hr = self.human.rhr
        for n in range(len(self.avg_hr)):
            t, hr, met = self.avg_hr[n]
            if hr != prev_hr:

                end_t = t + 2
                if len(self.avg_hr) - 1 > n:
                    # check the next one
                    next_t = self.avg_hr[n + 1][0]
                    if next_t < end_t:
                        end_t = next_t

                # Add transition steps
                t_steps = np.arange(t, end_t + t_step, t_step)
                hr_steps = np.linspace(prev_hr, hr, num=len(t_steps))
                temp_events.extend([(ts, hr_steps[i], met) for i, ts in enumerate(t_steps)])

            prev_hr = hr

        # Write the HR data to the output file
        self.avg_hr = temp_events
        if output_file_path is not None:
            with open(output_file_path, 'w') as hr_output_file:
                for t, hr, met in self.avg_hr:
                    hr_output_file.write("{},{},{}\n".format(t, hr, met))

        # If the visual flag is set plot the results.
        if self.visual:
            data = np.array(self.avg_hr)
            plt.figure()
            plt.plot(data[:, 0], data[:, 1])
            plt.xlabel("Time (min)")
            plt.ylabel("Average HR")
            plt.grid(True)
            plt.show()

        return 0  # return status


class ECGSimulation(object):
    """Simulate realistic ECG signals based on inputs."""

    def __init__(self, hr_list, visual=False):
        self.hr_vector = hr_list
        self.ecg_state = None
        self.x0 = np.array([1, 0, 0.04])
        self.visual = visual

    def synthetic_ecg(self, hr_mean=60.0, no_amp=0.0, start_t=0, stop_t=10):
        """
        :param hr_mean: Mean heart rate
        :param no_amp: Noise amplitude
        :param start_t: Signal start time in seconds
        :param stop_t: signal stop time in seconds
        :return: ECG signal array and corresponding time (sec) array
        """

        # Settings ----------------------------------------------------------- #
        Fs = 100.0                                 # sampling frequency (samples/sec)
        sfecg = Fs                                 # ECG sampling frequency [ Hertz]
        sfint = Fs                                 # ECG sampling frequency [ Hertz]
        hr_std = 1.0                               # Standard deviation of heart rate [1 beat per minute]
        lfhfratio = 0.5                            # LF/HF ratio [0.5]
        ti = np.radians([-60, -15, 0, 15, 90])     # P  Q  R  S  T   ti = angles of extrema degrees [radians]
        ai = np.array([1.2, -5, 30, -7.5, 0.75])   # ai = z-position of extrema [1.2 -5 30 -7.5 0.75]
        bi = np.array([0.25, 0.1, 0.1, 0.1, 0.4])  # bi = Gaussian width of peaks [0.25 0.1 0.1 0.1 0.4]
        min_rand(9843)                             # Seed the RNG
        # -------------------------------------------------------------------- #

        n = int(stop_t - start_t)  # time in sec

        # Adjust extrema parameters for mean heart rate
        hrfact = np.sqrt(hr_mean/60.0)
        hrfact2 = np.sqrt(hrfact)
        bi = hrfact*bi
        ti = np.array([hrfact2, hrfact, 1., hrfact, hrfact2])*ti
        q = np.round(sfint/sfecg)

        # frequency parameters for rr process flo and fhi are the Mayer waves
        # and respiratory rate respectively
        flo = 0.1
        fhi = 0.25
        flo_std = 0.01
        fhi_std = 0.01

        # Compute the RR-interval which is the time between successive R-peaks,
        # the inverse of this time interval gives the instantaneous heart rate.
        sfrr = 1.0  # sampling frequency
        w1 = 2*np.pi*flo
        w2 = 2*np.pi*fhi
        c1 = 2*np.pi*flo_std
        c2 = 2*np.pi*fhi_std
        sig2 = 1.
        sig1 = lfhfratio
        rr_mean = 60./hr_mean
        rr_std = 60.*hr_std/(hr_mean*hr_mean)

        df = sfrr/n
        w = np.arange(n)*2*np.pi*df
        dw1 = w - w1
        dw2 = w - w2

        hw1 = sig1*np.exp(-0.5*(dw1/c1)**2)/np.sqrt(2*np.pi*c1**2)
        hw2 = sig2*np.exp(-0.5*(dw2/c2)**2)/np.sqrt(2*np.pi*c2**2)
        hw = hw1 + hw2
        sw = (sfrr/2.)*np.sqrt(hw)

        if n % 2 == 0:
            ph0 = 2*np.pi*np.array([min_rand() for _ in range(int(n/2 - 1))])
            ph = np.hstack((0, ph0, 0, -np.flipud(ph0)))
        else:
            ph0 = 2*np.pi*np.array([min_rand() for _ in range(int(n/2))])
            ph = np.hstack((0, ph0, 0, -np.flipud(ph0)))[0:-1]

        swc = sw*np.exp(1j*ph)
        x = (1./n)*np.real(np.fft.ifft(swc))
        xstd = np.std(x)
        ratio = rr_std/xstd
        rr = rr_mean + x*ratio

        # Up-sample rr time series from 1 Hz to sfint Hz
        rr = np.interp(np.linspace(0, len(rr) - 1, len(rr)*sfint), range(len(rr)), rr)

        # make the rrn time series
        dt = 1/sfint

        # Solve the ODE system  shown below
        t_span = np.arange(start=dt, stop=len(x), step=dt)
        sol = odeint(derivsecgsyn, self.x0, t_span, args=(rr, sfint, ti, ai, bi))

        # Access the z vector of the solution (noiseless ECG amplitude)
        # (Down sample to required ecg_sf)
        x = sol[::q, :][:, 0]
        y = sol[::q, :][:, 1]
        z = sol[::q, :][:, 2]
        #self.x0 = np.array([x, y, z])

        # Scale the ECG signal to lie between ECG_AMP_MIN and ECG_AMP_MAX
        z_min = z.min(0)
        z_max = z.max(0)
        z_range = z_max - z_min
        z = (z - z_min)*(ECG_AMP_MAX - ECG_AMP_MIN)/z_range + ECG_AMP_MIN

        # Include additive uniformly distributed measurement noise
        noise = no_amp*(2*np.array([min_rand() for _ in range(len(z))]) - 1)
        syn_ecg = z + noise

        # Plot ECG vs Time (debugging)
        if self.visual:
            plt.figure()
            plt.plot(t_span + start_t, syn_ecg)
            plt.xlabel("Time [sec]")
            plt.ylabel("ECG Amplitude [mV]")
            plt.grid(True)
            plt.show()

        return syn_ecg, t_span

    def run_simulation(self, file_path=None):
        """Create a synthetic ECG signal for the given HR vector."""

        ecg_signal = []
        for n, hr_obj in enumerate(self.hr_vector):
            t, hr, met = hr_obj
            if n < len(self.hr_vector) - 1:
                t_start_sec = t*60
                t_stop_sec = self.hr_vector[n + 1][0]*60
                ecg_signal.append(self.synthetic_ecg(hr, ECG_NOISE_AMP, t_start_sec, t_stop_sec))

        if file_path is not None:
            with open(file_path, 'w') as ecg_output_file:
                for it, trace in enumerate(ecg_signal):
                    ecg_vector, t_vector = trace
                    if it == 0:
                        ecg_output_file.write("{},{}\n".format(t_vector[0], t_vector[1] - t_vector[0]))
                    ecg_output_file.write(','.join(map(str, ecg_vector)) + "\n")

        return 0

# Functions ------------------------------------------------------------------ #


def min_rand(seed=0):
    """Minimal random number generator of Park and Miller. Returns a uniform
    random deviate between 0.0 and 1.0. Set or reset idum to any integer value
    (except the unlikely value MASK) to initialize the sequence; idum must not
    be altered between calls for successive deviates in a sequence.
    Ref. Numerical Recipes in C 2nd ed.

    :param seed: Set or reset seed to any integer value (except the unlikely
                 value MASK) to initialize the sequence; seed must not be
                 altered between calls for successive deviates in a sequence.
    """

    # Define constants
    ia = 16807
    im = 2147483647
    am = (1.0/im)
    iq = 127773
    ir = 2836
    mask = 123459876

    # Only allow the generator to be seeded once
    if "seed" not in min_rand.__dict__:
        # XORing with MASK allows use of zero and other simple bit patterns for
        # seed.
        min_rand.seed = seed ^ mask

    # Compute idum=(IA*idum) % IM without over-flows by Schrage's method.
    k = min_rand.seed/iq
    min_rand.seed = ia*(min_rand.seed - k*iq) - ir*k

    if min_rand.seed < 0:
        min_rand.seed += im
    ans = am*min_rand.seed  # Convert to a floating result.

    return ans


def read_input_file(file_path, verbose=False):
    """ Read the input file.

    Read the input file and return a human object and a queue of activities.

    :param file_path: Path to the input file.
    :param verbose: If true, print summary of file contents.
    :return: Tuple containing a human object and priority queue of activities.
    """

    # Read the input file (INI format)
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(file_path)

    # Read the human attributes
    human = Human(verbose=verbose)
    human.age = config.getfloat('human', 'age')
    human.gender = config.get('human', 'gender')
    human.mass = config.getfloat('human', 'mass')
    human.rhr = config.getfloat('human', 'resting_heart_rate')
    human.hr = human.rhr
    human.vo2max = config.getfloat('human', 'vo2max')

    if verbose:
        print(human.description())

    # Create the queue of activities
    activity_queue = Queue.PriorityQueue()
    for section in config.sections():
        if "activity" in section:

            # Create the activity object
            new_activity = Activity(uid=uuid4().int)
            new_activity.start = config.getfloat(section, 'start')
            new_activity.duration = config.getfloat(section, 'duration')
            new_activity.type = config.get(section, 'type')
            new_activity.met = config.getfloat(section, 'met')

            # Create the event wrapper
            new_event = Event(START_ACTIVITY, new_activity)

            # Queue the event (start_time, event)
            activity_queue.put((new_activity.start, new_event))

            # Print data to the console if the verbose flag is set
            if verbose:
                print(new_activity.description())

    return human, activity_queue


def derivsecgsyn(x, t, rr, sf_int, ti, ai, bi):
    """This file provides dxdt = F(t,x) taking input parameters:

    Order of extrema: [P Q R S T]

    :param x: Initial conditions of x, y, and z
    :param t: Time vector
    :param rr: Time between successive R-peaks (inverse is instantaneous HR)
    :param sf_int:  Internal sampling frequency [Hertz]
    :param ti: angles of extrema [radians]
    :param ai: z-position of extrema
    :param bi: Gaussian width of peaks
    :return: dx/dt of F(t,x)
    """

    x0, y0, z0 = x
    ta = np.arctan2(y0, x0)
    a0 = 1.0 - np.sqrt(y0**2 + z0**2)
    ip = np.floor(t*sf_int)
    w0 = 2*np.pi/rr[ip]  # Equation (4) in paper
    fresp = 0.25
    zbase = 0.00015*np.sin(2*np.pi*fresp*t)
    dti = np.fmod(ta - ti, 2*np.pi)
    dx1dt = a0*x0 - w0*y0
    dx2dt = a0*y0 + w0*x0
    dx3dt = -np.sum(ai*dti*np.exp(-0.5*(dti/bi)**2)) - 1.0*(z0 - zbase)
    return np.array([dx1dt, dx2dt, dx3dt])


# NOT ACTUALLY USED #
def vo2_max(age, gender):
    """Return average VO2 Max (METs) based on age and gender.

    Tabular values obtained from: "Exercise Standards, A statement for
    healthcare professionals from the American Heart Association."

    :param age: Age of subject in years
    :param gender: String value containing either "male" or "female"
    :return: Normal VO2max value in METs
    """

    if gender[0] == 'm':
        if age <= 39:
            return 12
        if age <= 49:
            return 11
        if age <= 59:
            return 10
        if age <= 69:
            return 9
        else:
            return 8
    else:
        if age <= 39:
            return 10
        if age <= 49:
            return 9
        else:
            return 8


# ---------------------------------------------------------------------------- #


def main():
    """Main simulation script. Modify inputs here."""

    # -------- User Inputs and Settings -------------------------------------- #

    # The content of the input file describes characteristics of the human and
    # their different activities.
    input_file_path = "input_short.txt"

    # Output file paths
    ecg_file_path = "output_ecg.txt"
    hr_file_path = "output_hr.txt"

    # If the verbose flag is set details are printed to the console.
    verbose = True

    # If the visual flag is set simulation outputs are plotted.
    # Note: Each plot window must be closed to alow the simulation to continue.
    visual = True

    # ------------------------------------------------------------------------ #

    # Read input file and create the human object and queue of activities.
    human, fel = read_input_file(input_file_path, verbose=verbose)

    # Create the HR simulation object
    hr_sim = HeartSimulation(human, fel, verbose=verbose, visual=visual)

    # Run the simulation
    hr_sim.run_simulation(hr_file_path)

    # Generate the ECG signals (store in file for separate processing)
    ecg_sim = ECGSimulation(hr_sim.avg_hr, visual=visual)
    ecg_sim.run_simulation(file_path=ecg_file_path)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
