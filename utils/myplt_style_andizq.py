import matplotlib

TINY_SIZE = 8
SMALL_SIZE = 10
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.rcParams['axes.linewidth'] = 3.3 #2.5  
matplotlib.rcParams['axes.titlepad'] = 20
matplotlib.rcParams['axes.labelpad'] = 12
matplotlib.rcParams['xtick.major.width']=2.7#2.3
matplotlib.rcParams['ytick.major.width']=2.7#2.3
matplotlib.rcParams['xtick.minor.width']=2.2#1.8
matplotlib.rcParams['ytick.minor.width']=2.2#1.8
matplotlib.rcParams['xtick.major.size']=5.5#5.0
matplotlib.rcParams['ytick.major.size']=5.5#5.0
matplotlib.rcParams['xtick.minor.size']=3.5#3.0
matplotlib.rcParams['ytick.minor.size']=3.5#3.0
    
matplotlib.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=MEDIUM_SIZE+3)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE+2)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=MEDIUM_SIZE-1)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=MEDIUM_SIZE-1)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=MEDIUM_SIZE-1)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE-1)  # fontsize of the figure title
