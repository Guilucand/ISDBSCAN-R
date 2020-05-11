//
// Created by andrea on 06/04/20.
//

#include <cstdio>
#include <cstdlib>
#include "structs.h"

IsdbscanResult isdbscan(std::vector<kd_point<double>> dataset, int k, int batch_size, bool stratif, bool approximate);

std::vector<kd_point<double>> read_file(const char *file_name) {

    FILE *file = fopen(file_name, "r");

    std::vector<kd_point<double>> dataset;

    size_t buf_len = 512 * 1024;
    char *line = (char*)malloc(512 * 1024);

    int len = 0;
    while ((len = getline(&line, &buf_len, file)) != -1) {
        line[len] = '\0';
        kd_point<double> pt;
        pt.index = dataset.size();
        double val;
        int processed = 0, read;
        while (processed < buf_len && sscanf(line + processed, "%lf%n", &val, &read) == 1) {
            pt.features.push_back(val);
            processed += read;
        }
        dataset.push_back(pt);
    }

    free(line);
    fclose(file);

    return dataset;
}

void write_output(std::vector<int> result, const char *path) {
    FILE *out = fopen(path, "w");

    for (int x : result) {
        fprintf(out, "%d\n", x);
    }
    fclose(out);
}

int main(int argc, const char **argv) {

    if (argc < 4) return 1;

    bool stratif = argc > 4;

    auto dataset = read_file(argv[1]);

    auto result = isdbscan(dataset, atoi(argv[3]), dataset.size(), stratif, false);
    write_output(result.clusters, argv[2]);
    if (stratif) {
        write_output(result.layer, argv[4]);
    }

}
