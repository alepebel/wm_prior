#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for the confirmation bias experiment (using predcision code)

## Change the response system
## Make sure that the information is not overwritten and save a logfile with the subject information
## Allow independent sessions for the same subject
## Add some feedback (score table) / check what
## Add instructions

version = 1.0

#cd('/home/node2/Experiments/PreDCN-master/predcision')
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


from random import random
from numpy import sin, pi, cos
from scipy import signal, stats
from numpy.linalg import norm
from pandas import DataFrame


from psychopy import visual, logging, core, event,  gui, data
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from psychopy.preferences import prefs
from psychopy import prefs

#import importlib
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
monitor_features['Hz'] =  'auto' #144 #60 this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often

   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']


st_prop = st.stim_config(ifi) #Loading stim properties

fixation_gauss = visual.PatchStim(win=win, mask='gauss', size=1, pos=[0,0], sf=0, color='black',units='cm')
fixation_square = visual.PatchStim(win=win, size=1, pos=[0,0], sf=0, color='black',units='cm')
cursor = visual.Circle(win,radius=0.7,pos=[0,0],fillColor='white',units='cm') #
#cursor = visual.PatchStim(win=win, size=1, pos=[0,0], sf=0, color='black',units='cm')
resp = visual.Circle(win,radius=0.7,pos=[0,0],lineColor='white',units='cm')

mouse=event.Mouse(win=win,visible=False)

date = datetime.datetime.now().strftime('%d%m%Y%H%M')
filename = subj_id + '_' + date + '.csv' # Name + '_' +  

cnames = ['subj','trial','block','delay','fixated',
           'R','T_Angle','stim_sd',
           'choice_x','choice_y','choiceAngle','choiceR',
           'RT','MT',
           'ts_b','ts_f','ts_p','ts_d','ts_r','ts_e',
           'm_pos_x','m_pos_y','m_clock']
if log_on:
    with open(filename, 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(cnames)
            
      

# Generating cue stimuli
def gkern(l=5, sig=1.):
    """\
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


def routines():
    exp.quit_handler(win)
    fixation.draw()
    win.update()
    x_m,y_m=mouse.getPos()  


## Experiment design
nblock_reps = 2
block_priors = np.random.permutation(np.array([0, 45, 225])) # lets make 0 to be random distribution (assigning kappa = 0 later)
block_priors = np.repeat(block_priors, nblock_reps) ; #np.tile(block_priors, nblocks)
block_prior_kappa = 1

expInfo['block_priors'] = [block_priors, block_prior_kappa] # block priors and kappa

 
main_exp  = {}
main_exp['nblocks']     = len(block_priors) # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['trial_reps']  = 50 #1 is equal to 1 minute (2 trials x 3 repetitions)



stimList = []
for delay in [500, 750]: # how many ms
      for stim_sd in ['L', 'H']:
            stimList.append({'delay':delay, 'stim_sd':stim_sd}) #, 'n_reps': n_reps




corr_lotery = []
mouse = event.Mouse(visible = True)
mouse.setPos([0,0])

win.mouseVisible = False 


for thisBlock in range(main_exp['nblocks']): # iterate over blocks
      
    #instr.block_start(win)
    # lets trigger the beggining of the experiment
    
    #instr.block_ID(win, block_type, True) # explicit information about the block type
    inst = visual.TextStim(win, pos = [0,8])
    inst.wrapWidth = 20
    inst.height = 1
      

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
        
    kap = 200000 # set an impossible value to start the loop
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
        xstim,ystim = exp.toCar(st_prop['radius'],angle)           
        
        t_delay = thisTrial['delay'] 
        if thisTrial['stim_sd']  == 'L':
            stim_ID = gabor_L
        elif thisTrial['stim_sd']  == 'H':  
            stim_ID = gabor_H
        
        stim_ID.setPos([xstim,ystim])
                   
        #if sst: p_port.write(b'52')  # send trigger block init
        #if sst: p_port.write(b'00')  # send trigger
        fixx=1
        ts_b = time()
        fixation = fixation_gauss
        routines()
        mouse.clickReset()
        
        x,y=mouse.getPos()  
        while not ( abs(x) < 1 and abs(y) < 1):
            exp.quit_handler(win)
            x,y=gpos=mouse.getPos()
            cursor.setPos(gpos, operation='', log=None)                                        
            cursor.draw() 
            routines()
        
        mouse.setPos([0,0])        
        fixation = fixation_square
        fixation.setColor('black')

        routines()

        trialClock.reset()
        ts_f = time()

       
        for frameN in range(int(st_prop['FP'])):
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_fp.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  

            routines()

      
        ts_p = time()
        for frameN in range(st_prop['CUE_frames']):
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
               
            stim_ID.draw()
            routines()
       

        trialClock.reset()
        ts_d=time()

        for frameN in range(int(delay/ifi)):
            routines()
      
        mouse.setPos((0,0))
        pos = mouse.getPos()
        resp.setFillColor("black")
        resp.setPos(pos)
        fixation = fixation_gauss
        choice_trial=0
        event.clearEvents()
        fixation.setColor("black")
        if delay == 0:
            stim.draw()
        else:
            pass
        routines()
        
        m_pos_x=[]; m_pos_y=[]
        m_clock=[]
        
   
        trialClock.reset()
        ts_r = time()

        #reaction time
        x,y = mouse.getPos()
        rad=0
        frameN = 0
        while ( abs(x) < 1 and abs(y) < 1): # and frameN in range(RESPONSE_MAX)):
            if delay == 0:
                stim_ID.draw()
            else:
                pass
            
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_resphase.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                 
            x,y=pos = mouse.getPos()
            resp.setPos(pos)
            resp.draw()
            routines()
            frameN += 1
            pass
            
            
        rtime = trialClock.getTime()                        # response time
        
        # movement time
        trialClock.reset()
        frameN = 0            
        #keep waiting for the release
        while (mouse.getPressed()[0]==0): # and frameN in range(RESPONSE_MAX)):
            x,y=pos=mouse.getPos()
            resp.setPos(pos)
            t=trialClock.getTime()
            resp.draw()
            m_pos_x.append(x)
            m_pos_y.append(y)
            m_clock.append(t)
            if frameN == 0:     
               if sst: win.callOnFlip(p_port.write, tg_respinit.encode()) # send trigger
            if frameN == 1:     
               if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins             
            routines()
            frameN += 1
            pass
        
        rep_pos=pos
            
        angle_T = angle
        r_T = st_prop['radius']
        ts_e = time()                                       # end time                   
        mtime = trialClock.getTime()                        # movement time 
        
        if mouse.getPressed()[0]==1:
           # if sst: p_port.write(b'31')  # send trigger
            choice_pos=rep_pos
            choice_ang= float("%.2f" % exp.getAngle(rep_pos))
            choice_r=float("%.2f" % (norm(rep_pos)))
            err_ang= angle_T- choice_ang
            err_r=float("%.2f" % (st_prop['radius']-(choice_r)))# -*- coding: utf-8 -*-
           # if sst: p_port.write(b'00')  # send trigger
        else:
            choice_pos=array([np.nan,np.nan])
            choice_ang=np.nan  
            choice_r=np.nan
            err_ang=np.nan
            err_r=np.nan
        
        if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
        win.flip()
        if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
        win.flip()
                           
        if err_ang>180:
            err_ang=err_ang-360
            
        prec_angle.append(err_ang)    
        trial_data = [subj_id, n_trial,thisBlock,delay,fixx,
                       r_T,angle_T,thisTrial['stim_sd'],
                       choice_pos[0],choice_pos[1],
                       choice_ang,choice_r,
                       rtime,mtime,
                       ts_b,ts_f,ts_p,ts_d,ts_r,ts_e,
                       m_pos_x,m_pos_y,m_clock]
        if log_on:
            with open(filename, 'a') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=';',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
                
                spamwriter.writerow(trial_data)

        # Feedback      
        exp.say_msg(np.round(err_ang),0.3,win)
        #elist[n_trial] = trial_data # append new row data to list
        thr_trials_var.append(trial_data)
        resp.setPos(pos)
        #resp.setFillColor('white')
        resp.draw()
    
    block['data'] = thr_trials_var
    main_exp['Exp_blocks'][thisBlock] =  block
    expInfo['main_exp'] =  main_exp

    toFile(resultspath, expInfo) #saving file to disk
            #trialClocktimes.append(trial_times)
        # Get datafiles in pandas format and attack to main Exp.variabl
        
        
                       
    prec_angle=np.array(prec_angle)                           
    prec_angle= prec_angle[~np.isnan(prec_angle)]                           
    abs_err = np.mean(abs(prec_angle))

   # if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # send trigger
    win.flip()




       
def say_msg(message,duration,win):
    msgClock=core.Clock()
    msgText=visual.TextStim(win=win, ori=0,
        text=message,
        pos=[0,0], height=30, #1.5,
        color='black',units="cm")
    win.flip()

main_exp['monitor_features'] = monitor_features
main_exp['stim_props'] = st_prop
#main_exp['lotery'] = corr_lotery



expInfo['main_exp'] =  main_exp

print('Overall, %i frames were dropped.' % win.nDroppedFrames)


toFile(resultspath, expInfo) #saving file to disk    
    

#stim_input.close()

#clean up
event.clearEvents()
exp.exit_task(win)


