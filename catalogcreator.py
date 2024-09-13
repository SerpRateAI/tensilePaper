

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

def load_peaks_dataframe():
    df = pd.read_csv(f'peaks{day_number}.csv')
    df['datetime'] = df.init_arrival_time.apply(dates.num2date)
    df['obs_dt'] = df.datetime.apply(obspy.UTCDateTime)
    df['phone_number'] = df.arrival_hydrophone.apply(lambda p: int(p[1]))
    # we return the phone index to the array index which starts from 0
    # this returns the phone index back to where it was generated from
    # the minimum phone that it can be detected on is h3 so we return it
    # to 0 by taking 3
    # df['phone_index'] = df.phone_number - 1
    df['phone_index'] = df.phone_number - 3
    return df

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
    starttime = starttime - 0.2
    # endtime = starttime + 0.5
    endtime = starttime + 0.2
    trimmed = wf.copy().trim(starttime=starttime, endtime=endtime)
    return trimmed

def calc_aics(wf):
    """
    Uses obspy aic_simple to pick the start time of an event

    Parameters
    --------------------
    event : obspy.Stream
        an obspy stream with traces inside. the expected data
        will be only for a single event, not the whole data set

    Returns
    --------------------
    aic_t : list
        the times per hydrophone for each aic picked event
    aics : list
        the raw aics calculated for each event
    """

    # calculates aic score
    aics = [trigger.aic_simple(tr.data) for tr in wf]

    # finds minimum and returns index for aic scores
    diffs = [np.diff(aic, n=1) for aic in aics]
    maxes = [np.argmax(diff) for diff in diffs]
    

    # uses minimum index to retrieve the timestamp
    aic_t = [wf[n].times('matplotlib')[i] for n, i in enumerate(maxes)]

    # return aic_t, aics
    return maxes, aic_t, aics

def make_event(id):
    e = peaks_df.loc[id]

    # initial picking method for hydrophones that calculates the AIC
    # then uses the finite difference of the AIC to determine the
    # location of the first arrivals on each hydrophone
    # these are then ordered
    wf = get_event_waveforms(wf=waveforms, starttime=e.obs_dt)
    # mxs, aict, aics = calc_aics(wf=wf)
    aics = [trigger.aic_simple(tr.data) for tr in wf]

    for n, aic in enumerate(aics):
        np.savetxt(fname=f'aics/{id}_{n}.txt', X=aic)
        

    aict = np.array([dates.num2date(t) for t in aict])

    dts = []
    aict_diff = np.array([(aict[e.phone_index] - t).total_seconds() for t in aict])
    hphone_arrival_order = np.argsort(aict_diff)
    dts = aict_diff[hphone_arrival_order]
    dt = dts[1]
    
    # aict = np.array(aict)
    # A_aict = aict[e.phone_index]
    # B_aict = aict[aict != A_aict]

    # t_A = dates.num2date(A_aict)
    # t_Bs = [dates.num2date(t) for t in B_aict]

    # dts = np.array([(t_A - t_B).total_seconds() for t_B in t_Bs])
    # abs_dts = np.argmin(np.abs(dts))
    # # dt = t_Bs[abs_dts]
    # dt = dts[abs_dts]
    
    # t_A = dates.num2date(aict[hphone_arrival_order[0]])
    # t_Bplus = dates.num2date(aict[hphone_arrival_order[1]])
    
    # hphone_arrival_order = np.argsort(aict)
    # picking_method = 0
    

    # if hphone_arrival_order[0] + 2 != e.phone_index:
    #     # if the order doesn't match the original starting hydrophone
    #     # retry now with a hanning window put on the data set        
    #     wf_h = wf.copy().taper(type='hann', max_percentage=0.5)
    #     mxs, aict, aics = calc_aics(wf_h)
    #     hphone_arrival_order = np.argsort(aict)
    #     picking_method = 1


    # if hphone_arrival_order[0] + 2 != e.phone_index:
    #     # if the order doesn't match after the hanning window is applied
    #     # retry now, without the hanning window, but a smaller event
    #     # window
    #     starttime = e.obs_dt - 0.1
    #     endtime = e.obs_dt + 0.4
    #     wf_t = wf.copy().trim(starttime=starttime, endtime=endtime)
    #     mxs, aict, aics = calc_aics(wf_t)
    #     hphone_arrival_order = np.argsort(aict)
    #     picking_method = 2

    
    # if hphone_arrival_order[0] + 2 != e.phone_index:
    #     # if the order doesn't match after trimming the event window
    #     # set the arrival order to match the initial arrival pick
    #     # and randomly choose the phone above or below to calculate
    #     # the depth
    #     phones = np.array([0, 1, 2, 3])

    #     first_hydrophone = e.phone_index
    #     second_hydrophone = e.phone_index + 1 if e.phone_index < 3 else e.phone_index - 1
    #     remaining_hphones = list(set(phones).difference([first_hydrophone, second_hydrophone]))
        
    #     hphone_arrival_order = [first_hydrophone
    #                             ,second_hydrophone
    #                             ,remaining_hphones[0] # third
    #                             ,remaining_hphones[1] # fourth
    #                            ]

    #     mxs, aict, aics = calc_aics(wf)
    #     aict = [aict[i] for i in hphone_arrival_order]
    #     picking_method = 3

    # calculate dt for depth
    # t_A = dates.num2date(aict[hphone_arrival_order[0]])
    # t_B = dates.num2date(aict[hphone_arrival_order[1]])

    # dt = (t_A - t_B).total_seconds()
    # dt_pick_method = 0

    # if the dt that is calculated is larger than the maximum
    # travel time between the two hydrophones set the dt to
    # the maximum
    # if dt > 0.04:
    #     dt = 0.04
    #     dt_pick_method = 1
    # elif dt < -0.04:
    #     dt = -0.04
    #     dt_pick_method = -1
    # else:
    #     pass

    # calc relative depth
    # sign = hphone_arrival_order[0] - hphone_arrival_order[1]
    # dz_phone_label = 'h' + str(np.min([hphone_arrival_order[0], hphone_arrival_order[1]]) + 2)
    # hydrophone_depth = hydrophone_metadata[dz_phone_label]['depth']
    
    hydrophone_depth = hydrophone_metadata[e.arrival_hydrophone]['depth']
    sign = hphone_arrival_order[0] - hphone_arrival_order[1]
    relative_depth = 35 - 0.5 * dt * velocity_model * sign
    # relative_depth = 35 - 0.5 * dt * velocity_model
    depth = hydrophone_depth + relative_depth

    max_amp = wf[hphone_arrival_order[0]].data.max()
    # origin_time = obspy.UTCDateTime(t_A) - (relative_depth/velocity_model)
    origin_time = obspy.UTCDateTime(aict[e.phone_index]) - (relative_depth/velocity_model)
    
    event = {
        'id':id
        # arrival times decided by AIC
        # ,'t3':aict[hphone_arrival_order[0]]
        # ,'t4':aict[hphone_arrival_order[1]]
        # ,'t5':aict[hphone_arrival_order[2]]
        # ,'t6':aict[hphone_arrival_order[3]]
        ,'t3':aict[0]
        ,'t4':aict[1]
        ,'t5':aict[2]
        ,'t6':aict[3]
        # ,'aic3':aics[0]
        # ,'aic4':aics[1]
        # ,'aic5':aics[2]
        # ,'aic6':aics[3]
        # ,'first_hydrophone':hphone_arrival_order[0] + 2
        ,'first_hydrophone':hphone_arrival_order[0]
        # ,'second_hydrophone':hphone_arrival_order[1] + 2
        ,'second_hydrophone':hphone_arrival_order[1]
        ,'relative_depth':relative_depth
        ,'depth':depth
        # ,'picking_method':picking_method
        ,'dt':dt
        # ,'dt_pick_method':dt_pick_method
        ,'max_amp':max_amp
        ,'origin_time':origin_time
        ,'arrival_hydrophone':e.arrival_hydrophone
        ,'init_arrival_time':e.init_arrival_time
    }
    return event
    

if __name__ == '__main__':
    
    args = sys.argv
    day_number = args[1]

    swarm_starttime = config.swarm_starttime[day_number]
    swarm_endtime = config.swarm_endtime[day_number]

    hydrophone_metadata = config.hydrophone_metadata[day_number]

    paths = useful_variables.make_hydrophone_data_paths(borehole='a', year=2019, julian_day=day_number)
    
    waveforms = load.import_corrected_data_for_single_day(paths=paths)
    
    # remove first two hydrophone data since its bad
    # waveforms = waveforms[2:]

    peaks_df = load_peaks_dataframe()

    events = []
    for id in tqdm(peaks_df.index, desc='detecting events'):
        e = make_event(id=id)
        events.append(e)

    df = pd.DataFrame(events)
    df.to_csv(str(day_number)+'precision.csv')
    
    print('final number of events:', df.shape)
    
    print('first event time:', df.t3.min())
    print('last event time:', df.t4.max())

