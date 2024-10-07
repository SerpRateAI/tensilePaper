
import numpy as np
import obspy
from hydrophone_data_processing import load, useful_variables, plotting, signal_processing
import scipy.signal as signal
import pandas as pd
import matplotlib.dates as dates
import obspy.signal.trigger as trigger
import sys
import config
from tqdm import tqdm

args = sys.argv
day_number = args[1]

swarm_starttime = config.swarm_starttime[day_number]
swarm_endtime = config.swarm_endtime[day_number]

hydrophone_metadata = config.hydrophone_metadata[day_number]

paths = useful_variables.make_hydrophone_data_paths(borehole='a', year=2019, julian_day=day_number)
waveforms = load.import_corrected_data_for_single_day(paths=paths)
waveforms = waveforms[2:]
waveforms.filter(type='highpass', freq=50, zerophase=False, corners=1)

hydrophone_metadata = config.hydrophone_metadata['141']
velocity_model = 1750


# def load_peaks_dataframe():
peaks_df = pd.read_csv(f'peaks{day_number}.csv')
peaks_df['datetime'] = peaks_df.init_arrival_time.apply(dates.num2date)
peaks_df['obs_dt'] = peaks_df.datetime.apply(obspy.UTCDateTime)
peaks_df['phone_number'] = peaks_df.arrival_hydrophone.apply(lambda p: int(p[1]))
# we return the phone index to the array index which starts from 0
# this returns the phone index back to where it was generated from
# the minimum phone that it can be detected on is h3 so we return it
# to 0 by taking 3
# df['phone_index'] = df.phone_number - 1
peaks_df['phone_index'] = peaks_df.phone_number - 3

def get_event_waveforms(wf, starttime):
    """
    Returns a 1 second long trimmed event for the starttime and endtime
    
    Applies hanning window to edges to smooth to zero.

    Parameters
    -----------------
    starttime : obspy.UTCDatetime
        starttime for the event

    Returns
    -----------------
    data : obspy.Stream
        trimmed waveforms for 1 second event window
    """
    st = starttime - 0.2
    # endtime = starttime + 0.5
    et = starttime + 0.2
    trimmed = wf.copy().trim(starttime=st, endtime=et)
    return trimmed

def calc_depth(t1, t2, h1, h2):
    hydrophones = {'h3':0, 'h4':1, 'h5':2, 'h6':3}
    n1 = hydrophones[h1]
    n2 = hydrophones[h2]
    
    dt = t1 - t2
    
    # sign = (n1 - n2)/np.abs(n1 - n2)
    # sign = (n2 - n1)
    sign = (n1 - n2)
    if n1 - n2 > 0:
        h1 = h2
    # relative_depth = 35 - 0.5 * dt * velocity_model * sign
    # relative_depth = 35 - 0.5 * dt * velocity_model * sign
    relative_depth = 35 - 0.5 * dt * velocity_model * sign
    hydrophone_depth = hydrophone_metadata[h1]['depth']
    depth = hydrophone_depth + relative_depth
    return relative_depth, depth

def make_event(id):
    hydrophones = {'h3':0, 'h4':1, 'h5':2, 'h6':3}
    
    e = peaks_df.loc[id]
    datetime = e.datetime
    obs_dt = e.obs_dt
    phone_number = e.phone_number
    arrival_hydrophone = e.arrival_hydrophone

    n1 = hydrophones[arrival_hydrophone]

    wf = get_event_waveforms(wf=waveforms, starttime=e.obs_dt)
    
    tmaxes = []
    datetimes = []
    for n, tr in enumerate(wf):
                
        aic = trigger.aic_simple(tr.data)
        
        aicdiff = np.diff(aic)
        
        diffmax = np.argmax(aicdiff)+1
        
        tmax = tr.times()[diffmax]
        tmaxes.append(tmax)
        
        datetime = tr.times('matplotlib')[diffmax]
        datetimes.append(datetime)    

    
    if arrival_hydrophone == 'h3':
        second_hydrophone = 'h4'
        n2 = hydrophones[second_hydrophone]
        
        relative_depth, depth = calc_depth(t1=tmaxes[n1], t2=tmaxes[n2], h1=arrival_hydrophone, h2=second_hydrophone)
                                           
    elif arrival_hydrophone == 'h4':
        second_i = 'h3'
        second_j = 'h5'
        
        ni = hydrophones[second_i]
        nj = hydrophones[second_j]

        ti = tmaxes[ni]
        tj = tmaxes[nj]
        t1 = tmaxes[n1]

        dti = np.abs(t1 - ti)
        dtj = np.abs(t1 - tj)

        second_phone = np.argmin([dti, dtj])
        if second_phone == 0:
            second_hydrophone = 'h3'
            n2 = hydrophones[second_hydrophone]
        
        elif second_phone == 1:
            second_hydrophone = 'h5'
            n2 = hydrophones[second_hydrophone]

        else:
            raise ValueError('h4 is not calculating an index...')

        relative_depth, depth = calc_depth(t1=tmaxes[n1], t2=tmaxes[n2], h1=arrival_hydrophone, h2=second_hydrophone)
        
        
    elif arrival_hydrophone == 'h5':
        second_i = 'h4'
        second_j = 'h6'
        
        ni = hydrophones[second_i]
        nj = hydrophones[second_j]

        ti = tmaxes[ni]
        tj = tmaxes[nj]
        t1 = tmaxes[n1]

        dti = np.abs(t1 - ti)
        dtj = np.abs(t1 - tj)

        second_phone = np.argmin([dti, dtj])
        if second_phone == 0:
            second_hydrophone = 'h4'
            n2 = hydrophones[second_hydrophone]
        
        elif second_phone == 1:
            second_hydrophone = 'h6'
            n2 = hydrophones[second_hydrophone]

        else:
            raise ValueError('h4 is not calculating an index...')

        relative_depth, depth = calc_depth(t1=tmaxes[n1], t2=tmaxes[n2], h1=arrival_hydrophone, h2=second_hydrophone)
        
    elif arrival_hydrophone == 'h6':
        second_hydrophone = 'h5'
        n2 = hydrophones[second_hydrophone]
        
        relative_depth, depth = calc_depth(t1=tmaxes[n1], t2=tmaxes[n2], h1=arrival_hydrophone, h2=second_hydrophone)
        
    else:
        raise ValueError(f'this is not a hydrophone: {arrival_hydrophone}')

    # e = peaks_df.loc[id]
    # datetime = e.datetime
    # obs_dt = e.obs_dt
    # phone_number = e.phone_number
    # arrival_hydrophone = e.arrival_hydrophone
    event = {'id':id
            ,'datetime':datetime
            ,'obs_dt':obs_dt
            ,'phone_number':phone_number
            ,'arrival_hydrophone':arrival_hydrophone
            ,'relative_depth':relative_depth
            ,'depth':depth
            ,'max_amp':wf[n1].data.max()
            ,'second_hydrophone':second_hydrophone
            ,'t1':tmaxes[n1]
            ,'t2':tmaxes[n2]
            ,'origin_time':obs_dt - relative_depth/velocity_model}

    return event

if __name__ == '__main__':
    events = []
    for id in tqdm(peaks_df.index, desc='Using arrival as first'):
        e = make_event(id=id)
        events.append(e)

    df = pd.DataFrame(events)
    df.to_csv(f'{day_number}everything.csv')
    print(f'final number of events: {df.shape}')







    
    