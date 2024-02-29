import matplotlib.pyplot as plt
import h5py as h5
import caesar as cs
import numpy as np

filepath = "/Backup00/ejones/m50n512/dust_shielding/Tdust_m50n512_097.hdf5"
fptr = h5.File(filepath)

dust_temperature = {"on":[], "off":[]}
for group in fptr.keys():
    for dataset in fptr[group].keys():
        data = fptr[group][dataset][:]
        key = group.split("_")[-1].strip()
        dust_temperature[key].extend(data.tolist())

fig, ax = plt.subplots()

if np.allclose(dust_temperature["on"], dust_temperature["off"]):
    print("DUST SHIELDING HAS DONE FUCK ALL")

for key in dust_temperature.keys():
    ax.hist(dust_temperature[key], log=True, density=True, histtype="step", label=f"dust shielding = {key}")
    
fig.legend()
fig.tight_layout()
fig.savefig("dust_temperatures.png")