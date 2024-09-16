
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
# return df

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

def make_event(id):
    # event = {}
    e = peaks_df.loc[id]

    wf = get_event_waveforms(wf=waveforms, starttime=e.obs_dt)
    
    tmaxes = []
    for n, tr in enumerate(wf):
        np.savetxt(f'wfs/{id}_h{n}.txt', X=tr.data)
        aic = trigger.aic_simple(tr.data)
        np.savetxt(f'aics/{id}_h{n}.txt', X=aic)
        aicdiff = np.diff(aic)
        diffmax = np.argmax(aicdiff)+1
        tmax = tr.times()[diffmax]
        tmaxes.append(tmax)

    depths = []
    dts = []
    keys = []
    hydrophones = ['h3', 'h4', 'h5', 'h6']
    for n1, t1 in enumerate(tmaxes):
        for n2, t2 in enumerate(tmaxes):
            
            dt = t1 - t2
            dts.append(dt)
            
            sign = (n1 - n2)/np.abs(n1 - n2)
            relative_depth = 35 - 0.5 * dt * velocity_model * sign

            hydrophone_depth = hydrophone_metadata[hydrophones[n1]]['depth']
            depth = hydrophone_depth + relative_depth
            depths.append(depth)
            
            key = hydrophones[n1]+'_'+hydrophones[n2]
            keys.append(key)

    
    event = {}
    for n, key in enumerate(keys):
        # print(key)
        event[key+'_dt'] = dts[n]
        event[key+'_depth'] = depths[n]
        event['id'] = id

    # print(event.keys())
    return event
            
if __name__=='__main__':

    # peaks_df = load_peaks_dataframe()

    events = []
    for id in tqdm(peaks_df.index, desc='calc everything'):
        e = make_event(id=id)
        e['datetime'] = peaks_df.datetime.loc[id]
        # e['origin_time'] = e.obs_dt - (e.relative_depth / velocity_model) 
        events.append(e)

    df = pd.DataFrame(events)
    df.to_csv(str(day_number)+'everything.csv')
    print('final number of events:', df.shape)
