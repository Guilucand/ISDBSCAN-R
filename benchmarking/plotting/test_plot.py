import os
import blank_plot

name = "tsneLiver"  # "can473"
# name = "pcaLiver"  # "can473"
chname = "tsneLiver"  # "can473"
k = 45

input_name = f"../data/{name}.txt"
cluster_name = f"../data/{name}.{k}.clustering.txt"
chart_name = f"../data/{chname}.txt"

# Lazy clustering
if not os.path.exists(cluster_name):
    os.system(f"../../cmake-build-release/ISDBSCAN_Testing {input_name} {cluster_name} {k}")

blank_plot.show_plot(chart_name, cluster_name)
