#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for wm task. This is use to test healthy participants. 
# There are some catch events to maintain eyes on fixation point during trial


#pip install -U pip setuptools
#pip install tobii_research

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
tobii = False

# Initiating eyetracker
if tobii:
    import tobii_research as tr
    #import time

    found_eyetrackers = tr.find_all_eyetrackers()
    my_eyetracker = found_eyetrackers[0]
    print("Address: " + my_eyetracker.address)
    print("Model: " + my_eyetracker.model)
    print("Name (It's OK if this is empty): " + my_eyetracker.device_name)
    print("Serial number: " + my_eyetracker.serial_number)

''' Variable names
columns = ['device_time_stamp', 'system_time_stamp', 'left_gaze_direction_unit_vector', 
           'left_gaze_direction_validity', 'left_gaze_origin_position_in_hmd_coordinates', 
           'left_gaze_origin_validity', 'left_pupil_diameter', 'left_pupil_validity', 
           'left_pupil_position_in_tracking_area', 'left_pupil_position_validity', 
           'right_gaze_direction_unit_vector', 'right_gaze_direction_validity', 
           'right_gaze_origin_position_in_hmd_coordinates', 'right_gaze_origin_validity',
           'right_pupil_diameter', 'right_pupil_validity', 'right_pupil_position_in_tracking_area', 
           'right_pupil_position_validity']
'''


global_gaze_data = {}
global_gaze_data['left_pos'] = np.array([np.nan,np.nan])
global_gaze_data['right_pos'] = np.array([np.nan,np.nan])
global_gaze_data['left_pupil'] = np.nan
global_gaze_data['right_pupil'] = np.nan
global_gaze_data['left_pupil_validity'] = np.nan
global_gaze_data['right_pupil_validity'] = np.nan
global_gaze_data['device_time'] = np.nan
global_gaze_data['system_time'] = np.nan

if tobii:
    
    def gaze_data_callback(gaze_data): # this is a callback function that will return values
        global global_gaze_data
        global_gaze_data['left_pos'] = gaze_data['left_gaze_point_on_display_area']
        global_gaze_data['right_pos'] = gaze_data['right_gaze_point_on_display_area']
        global_gaze_data['left_pupil'] = gaze_data['left_pupil_diameter']
        global_gaze_data['right_pupil'] = gaze_data['right_pupil_diameter']
        global_gaze_data['left_pupil_validity'] = gaze_data['left_pupil_validity']
        global_gaze_data['right_pupil_validity'] = gaze_data['right_pupil_validity']
        global_gaze_data['device_time'] = gaze_data['device_time_stamp']
        global_gaze_data['system_time'] = gaze_data['system_time_stamp']


    # function to combine eye and left eye data
    def get_combined_eyes(gdata):
        combined_eyes = {}
        LPos = np.array(gdata['left_pos'])
        RPos = np.array(gdata['right_pos'])
        combined_eyes['EyesPos'] = np.nanmean([LPos,RPos], axis = 0)

        LPup = np.array(gdata['left_pupil'])
        RPup = np.array(gdata['right_pupil'])
        combined_eyes['EyesPup'] = np.nanmean([LPup,RPup], axis = 0)  
        return combined_eyes

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

mon, expInfo['monitor'] = exp.define_monitor(monitores[1]) # select the correct monitor

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
ring = visual.Circle(win,radius=10,pos=[0,0],edges = 60, lineColor='black',units='cm', fillColor = [0,0,0], interpolate = True)



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
cnames = ['subj','trial','block', 'prior','kappa','ct','delay','fix','keypressed',
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
            
      

# Generating cue stimuli with gaussian envelop
def gkern(l=5, sig=1.):
    ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(sig))
    kernel = np.outer(gauss, gauss)
    return kernel / np.sum(kernel)


gaussian_L = gkern(1000, 4)
gabor_L = visual.ImageStim(win,image=100*gaussian_L,size=(200, 200), pos=(0, 0), interpolate=True)

#gaussian_H = 2.23 * gkern(1000, 6) # 2.23 is the scaling factor to match the SD of the high contrast gabor (np.max(gaussianH)/np.max(gaussianL))
#gabor_H = visual.ImageStim(win,image=100*gaussian_H,size=(200, 200), pos=(0, 0), interpolate=True)

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
def circdist(angles1,angles2):       
    return np.angle(np.exp(1j*angles1)/np.exp(1j*angles2))

#block_priors = np.random.permutation(np.array([0, 45, 225])) # lets make 0 to be random distribution (assigning kappa = 0 later)
#block_priors = np.repeat(block_priors, nblock_reps) ; #np.tile(block_priors, nblocks)
nblock_reps = 1
prior_offset = 135 # difference in degrees between both angles
prior_1 = np.random.uniform(low=1, high=360, size=1)[0] # random initial position
prior_2 = (prior_1 + (np.random.choice([-1,1] , 1)[0]*prior_offset)) % 360  # summing in circular space (but randomly flipping whether the difference is positive or negative)
            #resp_init = exp.toCar(st_prop['radius'],np.rad2deg(z))
if prior_1 > 180: prior_1 = prior_1 -360
if prior_2 > 180: prior_2 = prior_2 -360
prior_1 = np.round(prior_1, 1)
prior_2 = np.round(prior_2, 1)

block_priors = [0, 0, prior_1, prior_1, prior_1, 0, 0, prior_2, prior_2, prior_2]
if subj_id == 'test':   block_priors = [0] #[0] # for test demo only one block
prior_kappa = 1 # kappa for the von mises distribution


expInfo['block_priors'] = [block_priors, prior_kappa] # block priors and kappa

 
main_exp  = {}
main_exp['nblocks']     = len(block_priors) # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['trial_reps']  = 25 # 25 * 2 delays = 50 
if subj_id == 'test':
    main_exp['trial_reps']  = 8 # 8 * 2 delays
    sst = False


# Initiating parallel port
if sst: 
    p_port = serial.Serial('COM3', 115200, timeout = 0 )
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')




stimList = []
for delay in [1000, 3000]: # how many ms 1000 and 3000
      for stim_sd in ['L']:
            stimList.append({'delay':delay, 'stim_sd':stim_sd}) #, 'n_reps': n_reps


mouse = event.Mouse(visible = True)
mouse.setPos([0,0])

win.mouseVisible = False 



if subj_id == 'test':
    instr.main_instructions(win)

if tobii:
    my_eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA, gaze_data_callback, as_dictionary=True)
    #my_eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)


for thisBlock in range(main_exp['nblocks']): # iterate over blocks
    print(thisBlock)
    if subj_id != 'test':  
        instr.block_start(win, thisBlock )
    # lets trigger the beggining of the experiment
    if sst: win.callOnFlip(p_port.write, tg_block.encode())
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

    block_prior_kappa = prior_kappa # kappa for the von mises distribution
    if block_priors[thisBlock] == 0: # if mean is 0, then we will generate a uniform distribution
        block_prior_kappa = 0 +1e-8 # adding value because perfect 0 gives an error in vonmisses function
        

    prec_angle=[]
    elist = [None] * nTrials
  
    thr_trials_var = [cnames]

    # save kappa of trial
    # save if there was correct fixation or not
    
    for n_trial, thisTrial in enumerate(trials): # iterate over trials
        trep = True
        # generating angle from vonmisses
        while trep:
            angle = vonmises.rvs(block_prior_kappa, 0, size=1)[0]
            angle = np.rad2deg(angle)
            angle = angle + block_priors[thisBlock] # shifting and correcting angles to be in the range of -180 to 180
            if angle > 180: angle = angle -360
            if angle < -180: angle = angle +360
            
            # mapping orientation in x y coordinates
            xstim,ystim = exp.toCar(st_prop['radius']+1.5,angle)  # 1.5 corrects for the different in radius between the cue and the ring
            
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
            #ct = False
            
            # decidign when the CT should appear
            ct_frames_interval =  [(ct_frames * 2), delay_frames-(ct_frames * 2)] # cant be smaller than ct duration, and leave some time from onset and offset of delay
            delay_ct = np.random.randint(low = ct_frames_interval[0], high= ct_frames_interval[1]) 
        
            
            EyesPos = np.array([0.0,0.0]) # initialising eye position    
            eye_lim = 0.65

            if tobii:
                eyepos = []
                while EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                # gdata = global_gaze_data
                    Eyes =  get_combined_eyes(global_gaze_data)
                    eyepos.append(Eyes['EyesPos']) #
                    if len(eyepos) > 500:
                        lastpos = eyepos[-30:]
                        EyesPos = np.nanmean(lastpos, axis=0)
                        #print(EyesPos)
                        if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                
            kfix = True # is the fixation correct?
                    
            if thisTrial['stim_sd']  == 'L':
                stim_ID = gabor_L
            elif thisTrial['stim_sd']  == 'H':  
                stim_ID = gabor_H
            
            stim_ID.setPos([xstim,ystim])
                    
            #if sst: p_port.write(b'52')  # send trigger block init
            #if sst: p_port.write(b'00')  # send trigger

            ts_b = time()
            routines(False)
            mouse.clickReset()
            trialClock.reset()
            ts_f = time()
            kfix = True
            if kfix == True:
                eyepos = []
                for frameN in range(int(st_prop['FP'])):
                    #print(frameN)
                    if frameN == 0:     
                        if sst: win.callOnFlip(p_port.write, tg_fp.encode()) # send trigger
                    if frameN == 1:     
                        if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                    if tobii:
                        Eyes = get_combined_eyes(global_gaze_data)
                        eyepos.append(Eyes['EyesPos'])
                        if frameN == int(st_prop['FP'])-1:
                            EyesPos = np.nanmean(eyepos, axis=0)
                            if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                            if  EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                                kfix = False
                                break       
                    routines(False)
    
            
                # gdata = global_gaze_data
                  #  Eyes =  get_combined_eyes(global_gaze_data)
                   # eyepos.append(Eyes['EyesPos']) #
                   # if len(eyepos) > 500:
                   #     lastpos = eyepos[-30:]
                   #     EyesPos = np.nanmean(lastpos, axis=0)
                   #     print(EyesPos)
                   #     if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                


            ts_p = time()

            if kfix == True:
                eyepos = []
                for frameN in range(st_prop['CUE_frames']):
                    if frameN == 0:     
                        if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
                    if frameN == 1:     
                        if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                    stim_ID.draw()
                    if tobii:
                        Eyes = get_combined_eyes(global_gaze_data)
                        eyepos.append(Eyes['EyesPos'])
                        if frameN == st_prop['CUE_frames']-1:
                            EyesPos = np.nanmean(eyepos, axis=0)
                            if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                            if  EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                                kfix = False
                                break       
                    routines(False)
        

            trialClock.reset()
            ts_d=time()

            if kfix == True:
                eyepos = []
                for frameN in range(delay_frames):
                    if tobii:
                        Eyes = get_combined_eyes(global_gaze_data)
                        eyepos.append(Eyes['EyesPos'])
                        if frameN == delay_frames-1: # in last frame, veriify whether the eyes are still in the fixation window
                            EyesPos = np.nanmean(eyepos, axis=0)
                            if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                            if  EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                                kfix = False
                                break       
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

            #keep waiting for the release
            # Draw the begining of the response phase

            z = np.random.uniform(low=0, high=359, size=1)[0] # random initial position
            resp_init = exp.toCar(st_prop['radius'],np.rad2deg(z))
            mouse.setPos(resp_init)

            rinit_time = trialClock.getTime()                        # response time
            kb.clock.reset()  # when you want to start the timer from
            kb.clearEvents()
            keypressed = None

            posresp = (0, 0) # initialising the position of the response

            if kfix == True:
                while (mouse.getPressed()[0]==0) and (keypressed == None ): # and frameN in range(RESPONSE_MAX)):
                    keys = kb.getKeys(keyList = ['space'], waitRelease=True, clear=True)
                    for thisKey in keys:
                        keypressed = thisKey.name
                        print(keypressed, thisKey.tDown, thisKey.rt)
                    x,y=pos=mouse.getPos()
                    posresp = exp.toCar(st_prop['radius'],np.rad2deg(x))     
                    resp.setPos(posresp)
                    t=trialClock.getTime()   
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
            #if kfix == True:
            resp_time = trialClock.getTime() - rinit_time                  # response time                       # response time
            rep_pos=posresp
                
            angle_T = angle
            r_T = st_prop['radius']
            ts_e = time()                                       # end time                   
            mtime = trialClock.getTime()                        # movement time 
            
            # calculating variables
            if (mouse.getPressed()[0]==1) and (kfix == True):
            # if sst: p_port.write(b'31')  # send trigger
                choice_pos=rep_pos
                choice_ang= float("%.2f" % exp.getAngle(rep_pos))
                choice_r=float("%.2f" % (norm(rep_pos)))
                err_ang = np.rad2deg(circdist( np.deg2rad(angle_T),np.deg2rad(choice_ang)))
                #err_ang= circdist( np.deg2rad(angle_T),np.deg2rad(choice_ang)) 
                err_r=float("%.2f" % (st_prop['radius']-(choice_r)))# -*- coding: utf-8 -*-
            # if sst: p_port.write(b'00')  # send trigger
            else:
                choice_pos=np.array([np.nan,np.nan])
                choice_ang=np.nan  
                choice_r=np.nan
                err_ang=np.nan
                err_r=np.nan

            
            
            if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
            ring.draw()
            routines(False)
            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
            ring.draw()
            routines(False)
                            
           # if err_ang>180:
           #    err_ang=err_ang-360
                
            #prec_angle.append(err_ang)    
            trial_data = [subj_id, n_trial,thisBlock,block_priors[thisBlock],prior_kappa, ct, t_delay,kfix,keypressed,
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
            if abs_err_angle >= 5 and abs_err_angle < 25:
                col_fb = 'orange'
            if abs_err_angle >= 25 and abs_err_angle < 60:
                col_fb = 'red'
            if abs_err_angle > 60:
                col_fb = 'black'

            for frameN in range(st_prop['FEEDBACK_frames']):
                if kfix == True:
                    ring.draw()
                    resp.draw()
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
                
                else:
                    msgText=visual.TextStim(win=win, ori=0,
                    text= 'MANTEN LOS OJOS EN EL PUNTO DE FIJACION',
                    pos=[0,-2], height=1,
                    color='red',units="cm").draw()

                routines(False)
            if kfix == True:
                trep = False # end trial and begin the next one

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
        
        
                       
    #prec_angle=np.array(prec_angle)                           
    #prec_angle= prec_angle[~np.isnan(prec_angle)]                           
    #abs_err = np.mean(abs(prec_angle))

   # if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # send trigger
    routines(False)

if tobii:
    my_eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)

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


