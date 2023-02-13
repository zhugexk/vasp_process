import re
import matplotlib.pyplot as plt


with open("C:/Users/11159/Desktop/TMP/Cu2Hg1Se4Sn1/BAND.dat", "r") as f:
    lines = f.readlines()

res = []
data = []

for line in lines:
    if len(re.findall(r"K-Path", line)) > 0:
        continue
    elif len(re.findall(r"NKPTS", line)) > 0:
        continue
    elif len(re.findall(r"Band-Index", line)) > 0:
        data = []
        res.append(data)
    elif len(re.findall(r"-*[0-9]+.[0-9]+", line)) > 0:
        data.append({"k_path": float(re.findall(r"-*[0-9]+.[0-9]+", line)[0]),
                     "energy_level": float(re.findall(r"-*[0-9]+.[0-9]+", line)[1])})

print(res)

with open("C:/Users/11159/Desktop/TMP/Cu2Hg1Se4Sn1/KLABELS", 'r') as f:
    lines = f.readlines()

k_labels = []
for line in lines:
    print(line)
    if len(re.findall(r"[0-9]+.[0-9]+", line)) > 0:
        print(re.findall(r"[A-Z|]+", line))
        label = re.findall(r"[A-Z|]+", line)[0]
        position = float(re.findall(r"[0-9]+.[0-9]+", line)[0])
        k_labels.append({"label": label, "position": position})

print(k_labels)


for data in res:
    plt.plot([d["k_path"] for d in data], [d["energy_level"] for d in data])
plt.grid(visible=True)

plt.xticks([k["position"] for k in k_labels], labels=[k["label"] for k in k_labels], ha='right')
plt.show()


