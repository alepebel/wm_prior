#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for wm task. This is use to test healthy participants. 
# There are some catch events to maintain eyes on fixation point during trial
3526#

version = 1.0

import numpy as np
from scipy.stats import vonmises
import matplotlib.pyplot as plt
import os, sys, ctypes, csv
import exp_func as exp
import stimuli as st
import instructions as instr
import serial
from time import sleep,time
import datetime
from psychopy.hardware import keyboard



from random import random
from numpy import sin, pi, cos
from scipy import signal, stats
from numpy.linalg import norm
from pandas import DataFrame


from psychopy import visual, logging, core, event,  gui, data
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from psychopy.preferences import prefs
from psychopy import prefs

import importlib
#importlib.reload(exp)

#prefs.general['winType'] =  'glfw'
#print(prefs)

Clock = core.Clock()

# Some general presets
event.globalKeys.clear() # implementing a global event to quit the program at any time by pressing ctrl+q
event.globalKeys.add(key='q', modifiers=['ctrl'], func=core.quit)


log_on=True


sst = False # Using parallel port to send triggers

if sst: 
    p_port = serial.Serial('/dev/tty.usbserial-BBTKUSBTTL', 115200, timeout = 0) # 'COM3', 115200, timeout = 0 this is for windows
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')


#Triggers

tg_fp =  '07'
tg_stim = '27'
tg_block = '61'
tg_resphase = '41'
tg_respinit = '35'
tg_resp = '31'
tg_zero = '00'

#prefs.hardware['audioLib']=['pyo'] # use Pyo audiolib for good temporal resolution
#from psychopy.sound import Sound # This should be placed after changing the audio library
# monitors.getAllMonitors()
#from general_settings import * # reads the variables in one script

# Collect subject data
expInfo, resultspath = exp.mainexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']

# Loading monitor definitions
monitores = st.monitor_def() 

#mon, expInfo['monitor'] = exp.define_monitor(monitores[1]) # select the correct monitor
mon, expInfo['monitor'] = exp.define_monitor(monitores[0]) # select the correct monitor

# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'deg' # units to define your stimuli
monitor_features['screen_id'] = 0 # 1 when using a extended display 
monitor_features['full']  = False
monitor_features['Hz'] =  60 #'auto' #144 #60 this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

kb = keyboard.Keyboard()
   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']


st_prop = st.stim_config(ifi) #Loading stim properties

fixation = visual.PatchStim(win=win, mask='gauss', size=1, pos=[0,0], sf=0, color=[-1,-1,-1],units='cm')
fixation_ct = visual.PatchStim(win=win, mask='gauss', size=1, pos=[0,0], sf=0, color=[-0.75,-0.75,-0.75],units='cm')


resp = visual.Circle(win,radius=0.7,pos=[0,0],fillColor='red',units='cm')
ring = visual.Circle(win,radius=10,pos=[0,0],edges = 60, lineColor='black',units='cm', interpolate = True)



mouse=event.Mouse(win=win,visible=False)

# output files paths (csv and psychopy pickle)
date = datetime.datetime.now().strftime('%d%m%Y%H%M')
if subj_id == 'test':
    filename = subj_id + '.csv'
    filepy = subj_id + '.psydat'
else:
    filename = subj_id + '_' + date + '.csv'
    filepy = subj_id + '.psydat'

filename = resultspath +'/'+ filename
filepy = resultspath +'/'+ filepy

# variables recorded
cnames = ['subj','trial','block', 'prior','ct','delay','fixated','keypressed',
           'R','T_Angle','stim_sd',
           'choice_x','choice_y','init_resp_angle','choiceAngle','choiceR',
           'RT','MT',
           'ts_b','ts_f','ts_p','ts_d','ts_r','ts_e',
           'm_pos_x','m_pos_y','m_clock']
if log_on:
    with open(filename, 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(cnames)
            
      
def mycirc_dist(a1, a2):
    a1 = np.deg2rad(a1); a2 = np.deg2rad(a2);
    diff_a = np.min([(a1 - a2) % (2*np.pi), (a2 - a1) % (2*np.pi)])
    return np.rad2deg(diff_a)

# Generating cue stimuli
def gkern(l=5, sig=1.):
    """
    creates gaussian kernel with side length `l` and a sigma of `sig`
    """
    ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(sig))
    kernel = np.outer(gauss, gauss)
    return kernel / np.sum(kernel)


gaussian_L = gkern(1000, 4)
gabor_L = visual.ImageStim(win,image=100*gaussian_L,size=(200, 200), pos=(0, 0), interpolate=True)

gaussian_H = 2.23 * gkern(1000, 6) # 2.23 is the scaling factor to match the SD of the high contrast gabor (np.max(gaussianH)/np.max(gaussianL))
gabor_H = visual.ImageStim(win,image=100*gaussian_H,size=(200, 200), pos=(0, 0), interpolate=True)

ct_frames = np.round(50/ifi) # duration of catch trial

def routines(dim):
    exp.quit_handler(win)      
    if dim: 
        fixation_ct.draw()
    else:
            fixation.draw()
    win.update()
    x_m,y_m=mouse.getPos()  

## Experiment design

#block_priors = np.random.permutation(np.array([0, 45, 225])) # lets make 0 to be random distribution (assigning kappa = 0 later)
#block_priors = np.repeat(block_priors, nblock_reps) ; #np.tile(block_priors, nblocks)
nblock_reps = 1
prior_offset = 135 # difference in degrees between both angles
prior_1 = np.random.uniform(low=0, high=359, size=1) # random initial position
prior_2 = (prior_1 + (np.random.choice([-1,1] , 1)*prior_offset)) % 360  # summing in circular space (but randomly flipping whether the difference is positive or negative)
            #resp_init = exp.toCar(st_prop['radius'],np.rad2deg(z))

#block_priors = np.random.permutation(np.array([0, 45, 225])) # lets make 0 to be random distribution (assigning kappa = 0 later)
#block_priors = np.repeat(block_priors, nblock_reps) ; #np.tile(block_priors, nblocks)

block_priors = [0, 0, prior_1, prior_1, prior_1, 0, 0, prior_2, prior_2, prior_2]
if subj_id == 'test':   block_priors = [0] # for test demo only one block
block_prior_kappa = 3

expInfo['block_priors'] = [block_priors, block_prior_kappa] # block priors and kappa

 
main_exp  = {}
main_exp['nblocks']     = len(block_priors) # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['trial_reps']  = 25 # 25 * 2 delays = 50 
if subj_id == 'test':
    main_exp['trial_reps']  = 25 # 8 * 2 delays # AMOUNT OF TRIALS


stimList = []
for delay in [1000, 3000]: # how many ms 1000 and 3000
      for stim_sd in ['L']:
            stimList.append({'delay':delay, 'stim_sd':stim_sd}) #, 'n_reps': n_reps


mouse = event.Mouse(visible = True)
mouse.setPos([0,0])

win.mouseVisible = False 

if subj_id == 'test':
    instr.main_instructions(win)


for thisBlock in range(main_exp['nblocks']): # iterate over blocks

    if subj_id != 'test':  
        instr.block_start(win, thisBlock )
    # lets trigger the beggining of the experiment
      

    if sst: win.callOnFlip(p_port.write, tg_mblock.encode())
    win.flip()
    if sst: win.callOnFlip(p_port.write, tg_zero.encode())
    win.flip()
    
    BlockClockStart = Clock.getTime() # block experiment time
    block = {} # dummy variable to save block data
    
    trialClock = core.Clock()

    trialClocktimes = np.array([]) # saving whole times here
    correct_seq = np.array([]) # saving seq. error responses per trial sequence

    # generate trials vector
    trials = data.TrialHandler(stimList, main_exp['trial_reps'], method='random')
    nTrials = trials.nTotal
    # generate orientations vector (constraining to the chosen parameters)

    if block_priors[thisBlock] == 0: # if mean is 0, then we will generate a uniform distribution
        block_prior_kappa = 0 +1e-8 # adding value because perfect 0 gives an error in vonmisses function
        
    kap = 200000 # set an impossible value to start the loop and veryfy that kappa is as expected
    while (kap < block_prior_kappa - 0.1) or (kap > block_prior_kappa + 0.1):
        angles = vonmises.rvs(block_prior_kappa, 0, size=nTrials)
        kap, loc, scl = vonmises.fit(angles,fscale=1)
    
    angles = np.rad2deg(angles)
    angles = angles + block_priors[thisBlock]

    prec_angle=[]
    elist = [None] * nTrials
  
    thr_trials_var = [cnames]

    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots(1, 1)
    #ax.hist((angles), density=True, bins='auto', histtype='stepfilled', alpha=0.2)
    #ax.legend(loc='best', frameon=False)
    #plt.show()


    for n_trial, thisTrial in enumerate(trials):  # will continue the staircase until it terminates!         
        # decision variable in this trial
        #instr.new_trial(win)
        angle = angles[n_trial]
        xstim,ystim = exp.toCar(st_prop['radius']+1.5,angle)           
        
        t_delay = thisTrial['delay'] 
        delay_frames = int(t_delay/ifi)

        # is this trial a catch trial?
        ct = False
        if subj_id == 'test':
            if np.random.randint(low = 0, high = 10) < 3: # for the demo increase the number of catch trials 
                ct = True
        else:
            if np.random.randint(low = 0, high = 10) == 0: 
                ct = True
        
        #ct = True # CHANGE TO TRUE LATER!!!!!!
        
        # decidign when the CT should appear
        ct_frames_interval =  [(ct_frames * 2), delay_frames-(ct_frames * 2)] # cant be smaller than ct duration, and leave some time from onset and offset of delay
        delay_ct = np.random.randint(low = ct_frames_interval[0], high= ct_frames_interval[1]) 


        if thisTrial['stim_sd']  == 'L':
            stim_ID = gabor_L
        elif thisTrial['stim_sd']  == 'H':  
            stim_ID = gabor_H
        
        stim_ID.setPos([xstim,ystim])
                   
        #if sst: p_port.write(b'52')  # send trigger block init
        #if sst: p_port.write(b'00')  # send trigger
        fixx=1
        ts_b = time()
        routines(False)
        mouse.clickReset()



        trialClock.reset()
        ts_f = time()

       
        for frameN in range(int(st_prop['FP'])):
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_fp.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  

            routines(False)
      
        ts_p = time()
        for frameN in range(st_prop['CUE_frames']):
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
               
            stim_ID.draw()
            routines(False)
       

        trialClock.reset()
        ts_d=time()

        
        for frameN in range(delay_frames):
            if ct and ((frameN >= delay_ct) and (frameN < delay_ct+ct_frames)):
                routines(True)
            else:
                routines(False)
      
        mouse.setPos((0,0))
        pos = mouse.getPos()
        choice_trial=0
        event.clearEvents()
        #fixation.setColor("black")
        
        m_pos_x=[]; m_pos_y=[]
        m_clock=[]     
   
        trialClock.reset()
        ts_r = time()
            
        # movement time
        trialClock.reset()

        #keep waiting for the release
        # Draw the begining of the response phase

        z = np.random.uniform(low=0, high=359, size=1) # random initial position
        resp_init = exp.toCar(st_prop['radius'],np.rad2deg(z))
        mouse.setPos(resp_init)

        rinit_time = trialClock.getTime()                        # response time
        kb.clock.reset()  # when you want to start the timer from
        kb.clearEvents()
        keypressed = None
        while (mouse.getPressed()[0]==0) and (keypressed == None ): # and frameN in range(RESPONSE_MAX)):
            keys = kb.getKeys(keyList = ['space'], waitRelease=True, clear=True)
            for thisKey in keys:
                keypressed = thisKey.name
                print(keypressed, thisKey.tDown, thisKey.rt)
            x,y=pos=mouse.getPos()
            posresp = exp.toCar(st_prop['radius'],np.rad2deg(x))     
            resp.setPos(posresp)
            t=trialClock.getTime()   
            #xstim,ystim = exp.toCar(st_prop['radius']+1.5,angle)           
            #stim_ID = gabor_L
            #stim_ID.setPos([xstim,ystim])
            #stim_ID.draw()       
            ring.draw()
            resp.draw()
            m_pos_x.append(x)
            m_pos_y.append(y)
            m_clock.append(t)
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_respinit.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins             
            routines(False)
            frameN += 1
            pass
        
        resp_time = trialClock.getTime() - rinit_time                  # response time                       # response time
        rep_pos=posresp
            
        angle_T = angle
        r_T = st_prop['radius']
        ts_e = time()                                       # end time                   
        mtime = trialClock.getTime()                        # movement time 
        
        if mouse.getPressed()[0]==1:
           # if sst: p_port.write(b'31')  # send trigger
            choice_pos=rep_pos
            choice_ang= float("%.2f" % exp.getAngle(rep_pos))
            choice_r=float("%.2f" % (norm(rep_pos)))
            err_ang= mycirc_dist(angle_T, choice_ang)
            err_r=float("%.2f" % (st_prop['radius']-(choice_r)))# -*- coding: utf-8 -*-
           # if sst: p_port.write(b'00')  # send trigger
        else:
            choice_pos=np.array([np.nan,np.nan])
            choice_ang=np.nan  
            choice_r=np.nan
            err_ang=np.nan
            err_r=np.nan

        
        
        if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
        routines(False)
        if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
        routines(False)


                           
        #if err_ang>180:
         #err_ang=err_ang-360
            
        prec_angle.append(err_ang)    
        trial_data = [subj_id, n_trial,thisBlock,block_priors[thisBlock], ct, t_delay,fixx,keypressed,
                       r_T,angle_T,thisTrial['stim_sd'],
                       choice_pos[0],choice_pos[1],z,
                       choice_ang,choice_r,
                       resp_time,mtime,
                       ts_b,ts_f,ts_p,ts_d,ts_r,ts_e,
                       m_pos_x,m_pos_y,m_clock]
        if log_on:
            with open(filename, 'a') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=';',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
                spamwriter.writerow(trial_data)

        # Feedback color scale

        abs_err_angle = np.abs(np.round(err_ang)) 
        if abs_err_angle < 5:
            col_fb = 'green'
        if abs_err_angle >= 5 and abs_err_angle < 20:
            col_fb = 'orange'
        if abs_err_angle >= 25 and abs_err_angle < 60:
            col_fb = 'red'
        if abs_err_angle > 60:
            col_fb = 'black'

        for frameN in range(st_prop['FEEDBACK_frames']):
            if ct == True:
                if keypressed == 'space':
                    msgTex=visual.TextStim(win=win, ori=0,
                    text= 'Detectado',
                    pos=[0,-2], height=1,
                    color='green',units="cm")
                    msgTex.draw()
                else:
                    msgText=visual.TextStim(win=win, ori=0,
                    text= 'No detectado',
                    pos=[0,-2], height=1,
                    color='red',units="cm").draw()
            if ct == False:
                if keypressed == 'space':
                    msgText=visual.TextStim(win=win, ori=0,
                    text= 'Error de deteccion',
                    pos=[0,-2], height=1,
                    color='green',units="cm").draw()
                else:
                    msgText=visual.TextStim(win=win, ori=0,
                    text= abs_err_angle,
                    pos=[0,-2], height=1.5,
                    color=col_fb,units="cm").draw()
            ring.draw()
            resp.draw()
            routines(False)

        #exp.say_msg(abs_err_angle,0.3, win)

        #elist[n_trial] = trial_data # append new row data to list
        thr_trials_var.append(trial_data)
        sleep(np.random.uniform(0.1,0.3)) # add some jitter time between trials

        #resp.setPos(pos)
        #resp.setFillColor('white')
        #resp.draw()
    
    block['data'] = thr_trials_var
    main_exp['Exp_blocks'][thisBlock] =  block
    expInfo['main_exp'] =  main_exp

    toFile(filepy, expInfo) #saving file to disk
            #trialClocktimes.append(trial_times)
        # Get datafiles in pandas format and attack to main Exp.variabl
        
        
                       
    prec_angle=np.array(prec_angle)                           
    prec_angle= prec_angle[~np.isnan(prec_angle)]                           
    abs_err = np.mean(abs(prec_angle))

   # if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # send trigger
    routines(False)


main_exp['monitor_features'] = monitor_features
main_exp['stim_props'] = st_prop
#main_exp['lotery'] = corr_lotery

expInfo['main_exp'] =  main_exp

print('Overall, %i frames were dropped.' % win.nDroppedFrames)

toFile(filepy, expInfo) #saving variables to disk    

instr.end_experiment(win)   
#clean up
event.clearEvents()
win.close()                                                              # eyetracker
core.quit
#exp.exit_task(win)



