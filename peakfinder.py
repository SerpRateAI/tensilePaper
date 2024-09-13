
import numpy as np
import obspy
from hydrophone_data_processing import load, useful_variables, plotting, signal_processing
import scipy.signal as signal
import pandas as pd
import matplotlib.dates as dates
import obspy.signal.trigger as trigger
import sys
import config
from Event import Event

hydrophone_metadata = config.hydrophone_metadata['141']

if __name__ == '__main__':
    import sys
    
    args = sys.argv
    day_number = args[1]
    peaks_loc = args[2]

    swarm_starttime = config.swarm_starttime[day_number]
    swarm_endtime = config.swarm_endtime[day_number]

    hydrophone_metadata = config.hydrophone_metadata[day_number]

    paths = useful_variables.make_hydrophone_data_paths(borehole='a', year=2019, julian_day=day_number)

    # loads data for all hydrophones
    # converts to pascals
    # flips the sign on hydrophone 3 if there it is borehole B due to wiring problem
    waveforms = load.import_corrected_data_for_single_day(paths=paths)
    print('loading data from:', paths)
    # filter and transform data

    ## trim data to be only for swarm times
    waveforms.trim(starttime=swarm_starttime, endtime=swarm_endtime)
    print('trimming data to be between:\n', swarm_starttime, '\n', swarm_endtime)
    
    ## 50hz high pass filter
    waveforms.filter(type='highpass', corners=1, zerophase=False, freq=50)
    print('applying 50Hz, corners=1, no zerophase, high pass filter')
    
    ## make copy for precision detection later
    waveforms_copy = waveforms.copy()
    print('making copy of data')

    ## square amplitude
    # omg this is where the magic is happening, this is why the max amp is so big its squared...
    for n, tr in enumerate(waveforms):
        waveforms[n].data = tr.data**2

    print('squaring amplitude of data')
    
    ## peak finder
    ## finds peaks in datat for each hydrophone
    peak_times = {}
    for n, tr in enumerate(waveforms):
        hydrophone_id = 'h' + str(n+1)

        # trim data for only selected data from hydrophone data
        tr_trim = tr.slice(starttime=hydrophone_metadata[hydrophone_id]['start']
                           ,endtime=hydrophone_metadata[hydrophone_id]['end'])
        t = tr_trim.times('matplotlib')
        print('trimming data for detection on', hydrophone_id, 'with starttime \n'
              ,hydrophone_metadata[hydrophone_id]['start']
              ,'\n and end time \n'
              ,hydrophone_metadata[hydrophone_id]['end']
             )

        # apply peak finding algorithm
        idx, props = signal.find_peaks(tr_trim.data, height=0.25, distance=250)
        print('finding peaks in squared data')

        # record initial event times for each peak detected for each hydrophone
        peak_times[hydrophone_id] = np.array(t[idx])

    ### 
    print(peak_times.keys())
    for k in peak_times.keys():
        print('hydrophone', k, 'number of events:', len(peak_times[k]))
    ## store initial picks in dataframe

    
    df_picks = pd.DataFrame()
    index_start = 0
    for k in peak_times.keys():
        init_arrivals = peak_times[k]
        n_events = init_arrivals.shape[0]
        index = np.arange(index_start, index_start+n_events, 1)
        rows = pd.DataFrame({'arrival_hydrophone':(k,)*n_events
                      ,'init_arrival_time':init_arrivals
                     }, index=index
                    )
        index_start = n_events
        print('hydrophone', k, 'number of events:', rows.shape)
        df_picks = pd.concat([df_picks, rows])
        
    print('storing initial peaks in dataframe')
    print('number of initial events detected:', df_picks.shape)
    print('first event:', dates.num2date(df_picks.init_arrival_time.min()))
    print('last event:', dates.num2date(df_picks.init_arrival_time.max()))
    print(f'writing peak dataframe to {peaks_loc}...')
    df_picks.to_csv(peaks_loc, index=False)