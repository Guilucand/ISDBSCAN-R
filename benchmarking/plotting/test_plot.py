import os
import blank_plot

for pc in range(5, 21, 5):
    for k in range(5, 50, 5):
        # name = "tsne2Liver"  # "can473"
        name = f"pca{pc}Liver"  # "can473"
        # chname = "tsne2Liver"  # "can473"
        chname = "umapLiver"  # "can473"

        input_name = f"../data/{name}.txt"
        cluster_name = f"../data/{name}.{k}.clusters.txt"
        layers_name = f"../data/{name}.{k}.layers.txt"
        chart_name = f"../data/{chname}.txt"

        use_stratif = True

        # Lazy clustering
        # if not os.path.exists(cluster_name):
        os.system(f"../../cmake-build-release/ISDBSCAN_Testing {input_name} {cluster_name} {k} "
                  f"{layers_name if use_stratif else ''}")

        blank_plot.show_plot(chart_name, cluster_name, f"results-{name}-{'strat' if use_stratif else 'nostrat'}", False)
        # while True:
        #     os.system("sleep 10")
        #     continue
