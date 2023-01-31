#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# A list with the main experimental functions
# Alexis Pérez-Bellido 01/2023

from psychopy import visual, logging, core, event,  gui, data
from psychopy import monitors
from psychopy import data, gui
from math import cos,sin, atan2,radians,degrees
from psychopy.tools.filetools import fromFile, toFile
from scipy.stats import norm
import os, sys
import numpy as np



def mainexp_subject_info(version): # loading staircase subject data   
    subjInfo = {'subj_id':'test'}
    # present a dialogue to change params
    
    sdataDlg = gui.DlgFromDict(subjInfo, title='Identifica participante')
    subj_id = subjInfo['subj_id']
    
    #fileName = subj_id
    #fileName = fileName+'.psydat' # the extension is required in order to be read by pyshcophy wrapper function "fromFile"
    # set path to results file 
    resultspath = os.path.join(os.getcwd(),'results') # os.sep
    filepy = os.path.join(resultspath,subj_id + '.psydat') # os.sep '
    try:  # try to get a previous parameters file           
        expInfo = fromFile(filepy)
        print("Participant" + subj_id + " data loaded")
        #expInfo = {'observer':subj_id, 'ExpVersion': version, 'guess':0.3, 'gender (M/F)': '?', 'age': 0, 'hand (L/R)': '?'}    
        # present a dialogue to change params    
        sdataDlg= gui.DlgFromDict(expInfo['subjInfo'], title='Orientation decisions: Threshold estimation', fixed=['ExpVersion'], order = [ 'ExpVersion', 'observer']) # fixed=['dateStr']    
        
        if not sdataDlg.OK:
            sys.exit('Experiment cancelled by the participant')
            core.quit
            
        expInfo['maindateStr'] = data.getDateStr() 
            
    except:  # if not there then use a default set
        
        myDlg = gui.Dlg(title="This participant does NOT exists!")
        myDlg.addText('This participant does NOT exists. Create new participant data (ok) or correct id (cancel)')
        ok_data = myDlg.show()
        
        if myDlg.OK:  
            
            print(ok_data)
            # lets create a subject from the scractch.
            expInfo = {}
            expInfo['maindateStr'] = data.getDateStr()
            subjInfo = {'observer':subj_id, 'ExpVersion': version, 'gender (M/F)': '?', 'age': 0, 'hand (L/R)': '?'}    
            # present a dialogue to change params            
            gui.DlgFromDict(subjInfo, title='Orientation decisions: Threshold estimation', fixed=['ExpVersion'], order = [ 'ExpVersion', 'observer']) # fixed=['dateStr']    
               
            fileName = subjInfo['observer']
            fileName = fileName+'.psydat' # the extension is required in order to be read by pyshcophy wrapper function "fromFile"
            expInfo['subjInfo'] = subjInfo        
            
        else:
            sys.exit('Experiment cancelled by the participant')
            
        
    #toFile(resultspath, expInfo) #saving file
    return expInfo, resultspath

    
def define_monitor(monitor):
    # Lets define the monitor characteristics
    screen = monitors.Monitor(name=monitor['monitor_name'] )
    screen.setSizePix(monitor['monitor_pixels'])
    screen.setWidth(monitor['monitor_width'] )
    screen.setDistance(monitor['distance2monitor'])
    screen.save()
    #monitors.getAllMonitors()
    print("The monitor " + screen.name  + " has "  + str(screen.getSizePix()) + " pixels, a width of " + str(monitor['monitor_width']) + " cm and subjects seats at " + str(monitor['distance2monitor']) + " cm"   )
    return screen, monitor



def create_window(monitor_features):
    # Lets open a window and define monitor used
    win = visual.Window(size = monitor_features['monitor'].getSizePix(), monitor = monitor_features['monitor'] , units=monitor_features['units'], screen=monitor_features['screen_id'],fullscr=monitor_features['full'])
    
    # win.waitBlanking = False # Only relevant for ubuntus equipment in order to prevent half of the frame rate detection
    # this makes imposible to record flip time
    
    Hz = monitor_features['Hz']
    real_ifi = win.getMsPerFrame(nFrames=60, showVisual=False, msg='', msDelay=0.0)
    if monitor_features['Hz'] == "auto": # if the monitor is not very reliable it might fail
        Hz = win.getActualFrameRate(nIdentical=20, nMaxFrames=80,nWarmUpFrames=10)
        Hz = np.round(Hz)
        ifi = real_ifi[0] # inter flip interval in order to calculate the time of each refresh
    
    else:
        ifi = 1000/monitor_features['Hz']
    # Threshold to evaluate dropped frames
    win.refreshThreshold = 1/Hz+ 0.004
    monitor_features['usedHz'] = Hz
    monitor_features['ifi'] = ifi
    monitor_features['real_ifi'] = real_ifi
    win.recordFrameIntervals = True
    win.mouseVisible = True
    return win, monitor_features


def draw_this_text(win, text, height):
    out = visual.TextStim(win, pos=[0, height], height = 0)
    out.wrapWidth = 20
   # out.setHeight = height
    out.text = text
    out.draw()

def exit_task(win):
    win.close()                                                              # eyetracker
    sys.exit(0)
    
def toCar(r,deg):
    x=r*cos(radians(deg))  #factor of 1cm==37.795 pixels
    y=r*sin(radians(deg))
    return (x,y)

def getAngle(v):
    a=np.arctan2(v[1],v[0])
    b=np.rad2deg(a)
    return b

#Function to say "x" message and report time
def say_msg(message,duration,win):
    msgClock=core.Clock()
    msgText=visual.TextStim(win=win, ori=0,
        text=message,
        pos=[0,0], height=1.5,
        color='black',units="cm")

    t=0
    msgClock.reset()
    while t<duration:
        t=msgClock.getTime()            
        msgText.draw()
        win.flip()
        
def instr_msg(message,win):
    msgText=visual.TextStim(win=win, ori=0,
    text=message,
    pos=[0,0], height=30, #1.5,
    color='black',units="pix")

    mouse.clickReset()
    while mouse.getPressed()[0]==0:  
        msgText.draw()
        win.flip()    

def len2(x):
    if type(x) is not type([]):
        if type(x) is not type(array([])):
            return -1
    return len(x)

def phase2(x):
    if not isnan(x):
        return phase(x)
    return nan

def quit_handler(win):
    if event.getKeys("escape"): exit_task(win)

def is_fix(fixated):
    return 1

 
# These two functions have been programmed by D. Linares.
def getResponse(win, keys, Clock):
    allKeys = event.getKeys(timeStamped=Clock)
    if len(allKeys) > 0:
        allK = [allKeys[0][0]]
        t = allKeys[0][1]
        for thisKey in allK:
            for k in keys:
                if thisKey == k:
                    return ([t, k])
                if thisKey in ["q", "escape"]:
                    print('cierre forzado')
                    win.close()
                    core.quit()

def waitResponse(win, keys):
    thisResp = None
    while thisResp == None:
        allKeys = event.waitKeys()
        for thisKey in allKeys:
            for k in keys:
                if thisKey == k:
                    return (k)
            if thisKey in ['q', 'escape']:
                print('cierre forzado')
                win.close()
                core.quit()
                
                

#subj = { 'mu' : 0, 'sd' : 0.3}
#xdat = np.linspace(-1, 1, num=90)
#ydat = norm.cdf(xdat, loc=subj['mu'], scale=subj['sd'])
# plot subject decision kernel

#plt.plot(xdat, ydat,linestyle='-')
                  
def simulated_subject(intensity, subj):
    ydat = norm.cdf(intensity, loc=subj['mu'], scale=subj['sd'])
    rand_v = np.random.uniform(low=0.0, high=1.0)
    dv = 1 if ydat >= rand_v else 0 # decision rule
    return dv # return decision variable

