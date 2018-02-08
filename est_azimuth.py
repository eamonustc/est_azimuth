#!/home/eamon/anaconda2/bin/python
import numpy as np
import obspy
from obspy.signal.cross_correlation import correlate
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'font.size': 14})

deg2rad = np.pi / 180.0
rad2deg = 180.0 / np.pi

# reading data
st_e = obspy.read('/home/eamon/myanmar/main_shock/meas_azimu/EM.M012.HLE.SAC')
st_n = obspy.read('/home/eamon/myanmar/main_shock/meas_azimu/EM.M012.HLN.SAC')
st_z = obspy.read('/home/eamon/myanmar/main_shock/meas_azimu/EM.M012.HLZ.SAC')
data_e = st_e[0].data
data_n = st_n[0].data
data_z = st_z[0].data
# add random noise to synthetics
scale = np.min([np.max(abs(data_e)), np.max(abs(data_n)), np.max(abs(data_z))])
time_e = np.arange(0, st_e[0].stats.npts, 1) * st_e[0].stats.delta
time_n = np.arange(0, st_n[0].stats.npts, 1) * st_n[0].stats.delta
time_z = np.arange(0, st_z[0].stats.npts, 1) * st_z[0].stats.delta
# tp and ts
#tp = st_z[0].stats.sac.t1
tp = st_z[0].stats.sac.t2 + 300.0 - 5.0
#ts = st_z[0].stats.sac.t2
ts = st_z[0].stats.sac.t2 + 300.0 + 12.0
print tp, ts
# estimate azimuth for each window
sampling_rate = 1.0 / st_e[0].stats.sampling_rate
win_len = 20  # window length for correlation
num_windows = int((ts-tp) / sampling_rate / win_len)
start_win = int(tp / sampling_rate)

azmu = []
coherence = []
t_azmu = []
for i_win in np.arange(0, num_windows, 1):
    corr_ez = np.mean(data_z[start_win+i_win*win_len : start_win+(i_win+1)*win_len] * data_e[start_win+i_win*win_len : start_win+(i_win+1)*win_len])
    corr_nz = np.mean(data_z[start_win+i_win*win_len : start_win+(i_win+1)*win_len] * data_n[start_win+i_win*win_len : start_win+(i_win+1)*win_len])
    corr_zz = np.mean(data_z[start_win+i_win*win_len : start_win+(i_win+1)*win_len] * data_z[start_win+i_win*win_len : start_win+(i_win+1)*win_len])
    one_azmu = np.arctan(abs(corr_ez) / abs(corr_nz))
    if (corr_nz>0) & (corr_ez>0):
        one_azmu = one_azmu
    elif (corr_nz<0) & (corr_ez>0):
        one_azmu = np.pi - one_azmu
    elif (corr_nz<0) & (corr_ez<0):
        one_azmu = np.pi + one_azmu
    else:
        one_azmu = 2.0*np.pi - one_azmu
    R = corr_zz / np.sqrt(corr_ez**2 + corr_nz**2)
    A = R * np.cos(one_azmu)
    B = R * np.sin(one_azmu)
    one_coherence = 1.0 - np.mean((data_z[start_win+i_win*win_len : start_win+(i_win+1)*win_len] - A*data_n[start_win+i_win*win_len : start_win+(i_win+1)*win_len] - B*data_e[start_win+i_win*win_len : start_win+(i_win+1)*win_len])**2) / corr_zz
    coherence.append(one_coherence)
    one_azmu = one_azmu * rad2deg
    azmu.append(one_azmu)
    t_azmu.append((start_win+i_win*win_len)*sampling_rate)

# convert to numpy array
coherence = np.array(coherence)
azmu = np.array(azmu)
t_azmu = np.array(t_azmu)

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8.5, 5))
ax1 = axes.flat[0]
#ax1.set_axis_off()
a1p1, = ax1.plot(time_e, data_e, 'k-', linewidth=0.8)
a1p2, = ax1.plot(time_n, data_n, 'b-', linewidth=0.8)
a1p3, = ax1.plot(time_z, data_z, 'r-', linewidth=0.8)
print tp, ts
ax1.axvline(x=tp, linestyle='--', linewidth=4.0, color='green')
ax1.axvline(x=ts, linestyle='--', linewidth=4.0, color='green')
ax1.legend([a1p1, a1p2, a1p3], ["E", "N", "Z"])
ax1.set_xlim([t_azmu[0], t_azmu[-1]])
ax1.set_ylim(-80000, 80000)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude (cm/s)')

ax2 = axes.flat[1]
#ax2.set_axis_off()
sc = ax2.scatter(t_azmu, azmu, c=coherence, s=40, cmap='PuRd', vmin=0, vmax=1)
ax2.axvline(x=tp, linestyle='--', linewidth=4.0, color='green')
ax2.axvline(x=ts, linestyle='--', linewidth=4.0, color='green')
plt.xlim([t_azmu[0], t_azmu[-1]])
plt.ylim([0, 360])
plt.xlabel('Time (s)')
plt.ylabel('Azimuth')

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.2)
cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cbar = fig.colorbar(sc, cax=cb_ax)
cbar.set_label('Predicted Coherence')
plt.show()
