from psychopy import visual, event, core
import random
import numpy as np
# Experiment instructions




def main_instructions_ieeg(win):
       
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.text = "Práctica: \
    En este experimento tendrás que recordar la posición de un estimulo que aparecera brevemente alrededor del punto de fijacion.\
    Después desaparacerá, y moviendo el ratón lateralmente habrás de ajustar la posición de una marca roja para que coincida con la posición recordada.\
     \
    Cuando estés mas o menos segura de la posición, harás clic con el ratón y verás por cuantos grados de ángulo te has equivocado"
    
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    inst.text = "ES MUY IMPORTANTE QUE TRATES DE MANTENER LOS OJOS EN EL PUNTO DE FIJACIÓN DURANTE CADA TRIAL YA QUE SI EL EYE-TRACKER DETECTA QUE NO ESTAS MIRANDO AL CENTRO EL TRIAL SE REINICIARÁ\
    Al final de cada bloque puedes descansar lo que necesites. También si estas cansada puedes mover los ojos un poco al final de cada trial. \
    El experimento tiene una duración de 6 bloques de 6 minutos aproximadamente cada uno. Tomatelo con calma y gracias por participar!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return



def main_instructions(win):
       
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.text = "Práctica: \
    En este experimento tendrás que recordar la posición de un estimulo que aparecera brevemente alrededor del punto de fijacion.\
    Después desaparacerá, y moviendo el ratón lateralmente habrás de ajustar la posición de una marca roja para que coincida con la posición recordada.\
     \
    Cuando estés mas o menos segura de la posición, harás clic con el ratón y verás por cuantos grados de ángulo te has equivocado"
    
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    inst.text = "Práctica: \
    Además de recordar la posición tendrás que estar atenta al punto de fijación porque invariablemente este parpaderará sutílmente.\
    SI DETECTAS QUE EL PUNTO DE FIJACIÓN HA PARPADEADO, en lugar de reportar la orientación del estímulo, EN LA FASE DE RESPUESTA HABRÁS DE PULSAR LA BARRA ESPACIADORA.\
    Cuando reportes tu respuesta, recibirás feedback sobre si te has equivocado"
    inst.height = 0.7
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Pulsa spacio para continuar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    
    inst.text = "Práctica: \
    ES MUY IMPORTANTE QUE TRATES DE MANTENER LOS OJOS EN EL PUNTO DE FIJACIÓN DURANTE CADA TRIAL. \
    Al final de cada bloque puedes descansar lo que necesites. También si estas cansada puedes mover los ojos un poco al final de cada trial. \
        \
    El experimento tiene una duración de 10 bloques de 6 minutos aproximadamente cada uno. Tomatelo con calma y gracias por participar!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return

def block_start(win, nblock,  lenblocks):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    nextt = visual.TextStim(win, pos = [0,-6])
    nextt.height = 0.7
    nextt.color = "black"
    
    inst.text = 'Comenzamos bloque ' + str(nblock+1) + ' / ' +  str(lenblocks)
        
    nextt.text = "Cuando estes listo/a, pulsa espacio para comenzar "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)
    return

    
def new_trial(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    inst.text = "Nuevo trial"
    inst.draw()
    win.flip()
    core.wait(0.75)
    return



    

def block_ID(win, rep, autodraw):
    inst = visual.TextStim(win, pos = [0,8])
    inst.wrapWidth = 20
    inst.height = 0.9
    inst.color = 'red'
    if rep == 'repeat':
        inst.text = "Secuencias repetidas"
    else:
        inst.text = "Secuencias parecidas"
    inst.autoDraw = autodraw   
    inst.draw()  
    
    
    


def end_experiment(win):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Final del experimento!! Avisa al investigador"     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return


def end_experiment_lot(win, corr_lotery, nblocks):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = 'Final de esta parte del experimento!! El participante ha ganado ' + str(sum(corr_lotery)) + ' puntos de ' + str(nblocks) + '. Avisa al investigador'     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return

    
    

def lotery(win, block, ifi):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Loteria de trials. Pulsa -espacio- para comenzar el sorteo."     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    allKeys = []
    number_text = visual.TextStim(win, pos=[0,4], height = 2)
    inst.text = "Loteria de trials. \
    Pulsa -espacio- para parar el sorteo."     
    n_trials = len(block['data'])
    
    time_number = round(100/ifi) # time in frames that each lotery number is displayed on the screen
    
    
    allKeys = []
    while allKeys != ['space']:
        #angle = np.random.uniform(0,1,1)
        trial_ix =  random.randint(0,n_trials-1)
        corr = block['data'].loc[trial_ix,'correct']
        col = 'red' if corr == -1 else 'green'
        number_text.text = trial_ix
        number_text.color = col
        
        for iframe in  range(time_number):
             inst.draw()
             number_text.draw()
             win.flip()
        
        allKeys = event.getKeys()
        
    allKeys = []
    
    inst.text = "Genial, has ganado! Pulsa -espacio- para continuar"
    if corr == -1:
        inst.text = "Oh! Qué mala suerte. Pulsa -espacio- para continuar"
        

    inst.draw()
    number_text.draw()
    win.flip()
        
    event.waitKeys(keyList = ["space"]) 
    return corr 
        
        
        