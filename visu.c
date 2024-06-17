#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct data {
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double *values;
};

int read_data(struct data *data, const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  if(!fp) {
    printf("Could not open binary data file '%s'\n", filename);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fread(&data->nx, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->ny, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->nz, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->xmin, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->xmax, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->ymin, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->ymax, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->zmin, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->zmax, sizeof(double), 1, fp) == 1);
  if(ok) {
    int N = data->nx * data->ny * data->nz;
    if(N <= 0) {
      printf("Invalid number of data points %d\n", N);
      ok = 0;
    }
    else {
      data->values = (double*)malloc(N * sizeof(double));
      if(!data->values) {
        printf("Could not allocate data (%d doubles)\n", N);
        ok = 0;
      }
      else {
        ok = (fread(data->values, sizeof(double), N, fp) == N);
      }
    }
  }
  fclose(fp);
  if(!ok) {
    printf("Error loading data from '%s'\n", filename);
    return 1;
  }
  return 0;
}

void free_data(struct data *data)
{
  free(data->values);
}

double get(const struct data *data, int m, int n, int p)
{
  return data->values[data->ny * data->nx * p + data->nx * n + m];
}

int write_data_gmsh(const struct data *data, const char *filename, int step)
{
  double delta_x = (data->nx > 1) ?
    ((data->xmax - data->xmin) / (data->nx - 1)) : 0;
  double delta_y = (data->ny > 1) ?
    ((data->ymax - data->ymin) / (data->ny - 1)) : 0;
  double delta_z = (data->nz > 1) ?
    ((data->zmax - data->zmin) / (data->nz - 1)) : 0;

  FILE *fp = fopen(filename, step == 0 ? "w" : "a");
  if(!fp) {
    printf("Could not open Gmsh data file '%s'\n", filename);
    return 1;
  }
  printf("Writing '%s' (step %d)...\n", filename, step);
  if(step == 0) {
    fprintf(fp, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    fprintf(fp, "$Nodes\n%d\n", data->nx * data->ny * data->nz);
    for(int p = 0; p < data->nz; p++) {
      for(int n = 0; n < data->ny; n++) {
        for(int m = 0; m < data->nx; m++) {
          int t = data->ny * data->nx * p + data->nx * n + m;
          double x = data->xmin + m * delta_x;
          double y = data->ymin + n * delta_y;
          double z = data->zmin + p * delta_z;
          fprintf(fp, "%d %g %g %g\n", t, x, y, z);
        }
      }
    }
    fprintf(fp, "$EndNodes\n$Elements\n");
    if(data->nx == 1) {
      fprintf(fp, "%d\n", (data->ny - 1) * (data->nz - 1));
      int e = 1;
      for(int p = 0; p < data->nz - 1; p++) {
        for(int n = 0; n < data->ny - 1; n++) {
          int t1 = data->ny * data->nx * p + data->nx * n;
          int t2 = data->ny * data->nx * (p + 1) + data->nx * n;
          int t3 = data->ny * data->nx * (p + 1) + data->nx * (n + 1);
          int t4 = data->ny * data->nx * p + data->nx * (n + 1);
          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    }
    else if(data->ny == 1) {
      fprintf(fp, "%d\n", (data->nx - 1) * (data->nz - 1));
      int e = 1;
      for(int p = 0; p < data->nz - 1; p++) {
        for(int m = 0; m < data->nx - 1; m++) {
          int t1 = data->ny * data->nx * p + m;
          int t2 = data->ny * data->nx * (p + 1) + m;
          int t3 = data->ny * data->nx * (p + 1) + (m + 1);
          int t4 = data->ny * data->nx * p + (m + 1);
          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    }
    else if(data->nz == 1) {
      fprintf(fp, "%d\n", (data->nx - 1) * (data->ny - 1));
      int e = 1;
      for(int n = 0; n < data->ny - 1; n++) {
        for(int m = 0; m < data->nx - 1; m++) {
          int t1 = data->nx * n + m;
          int t2 = data->nx * (n + 1) + m;
          int t3 = data->nx * (n + 1) + (m + 1);
          int t4 = data->nx * n + (m + 1);
          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    }
    else {
      fprintf(fp, "%d\n", (data->nx - 1) * (data->nz - 1) * (data->ny - 1));
      int e = 1;
      for(int p = 0; p < data->nz - 1; p++) {
        for(int n = 0; n < data->ny - 1; n++) {
          for(int m = 0; m < data->nx - 1; m++) {
            int t1 = data->ny * data->nx * p + data->nx * n + m;
            int t2 = data->ny * data->nx * p + data->nx * (n + 1) + m;
            int t3 = data->ny * data->nx * p + data->nx * (n + 1) + (m + 1);
            int t4 = data->ny * data->nx * p + data->nx * n + (m + 1);
            int t5 = data->ny * data->nx * (p + 1) + data->nx * n + m;
            int t6 = data->ny * data->nx * (p + 1) + data->nx * (n + 1) + m;
            int t7 = data->ny * data->nx * (p + 1) + data->nx * (n + 1) + (m + 1);
            int t8 = data->ny * data->nx * (p + 1) + data->nx * n + (m + 1);
            fprintf(fp, "%d 5 2 1 1 %d %d %d %d %d %d %d %d\n", e++,
                    t1, t2, t3, t4, t5, t6, t7, t8);
          }
        }
      }
    }
    fprintf(fp, "$EndElements\n");
  }
  fprintf(fp, "$NodeData\n1\n\"%s\"\n1\n%d\n3\n%d\n1\n%d\n",
          filename, step, step, data->nx * data->ny * data->nz);
  for(int p = 0; p < data->nz; p++) {
    for(int n = 0; n < data->ny; n++) {
      for(int m = 0; m < data->nx; m++) {
        int t = data->ny * data->nx * p + data->nx * n + m;
        fprintf(fp, "%d %g\n", t, get(data, m, n, p));
      }
    }
  }
  fprintf(fp, "$EndNodeData\n");
  fclose(fp);
  return 0;
}

int write_data_csv(const struct data *data, const char *filename, int step)
{
  FILE *fp = fopen(filename, step == 0 ? "w" : "a");
  if(!fp) {
    printf("Could not open CSV data file '%s'\n", filename);
    return 1;
  }
  printf("Writing '%s' (step %d)...\n", filename, step);
  for(int p = 0; p < data->nz; p++) {
    for(int n = 0; n < data->ny; n++) {
      for(int m = 0; m < data->nx; m++) {
        fprintf(fp, "%g;", get(data, m, n, p));
      }
    }
  }
  fprintf(fp, "\n");
  fclose(fp);
  return 0;
}

void closest_index(const struct data *data, double x, double y, double z,
                   int *mc, int *nc, int *pc)
{
  int m = (int)((x - data->xmin) / (data->xmax - data->xmin) * data->nx);
  int n = (int)((y - data->ymin) / (data->ymax - data->ymin) * data->ny);
  int p = (int)((z - data->zmin) / (data->zmax - data->zmin) * data->nz);
  if(m < 0) m = 0;
  else if(m > data->nx - 1) m = data->nx - 1;
  if(n < 0) n = 0;
  else if(n > data->ny - 1) n = data->ny - 1;
  if(p < 0) p = 0;
  else if(p > data->nz - 1) p = data->nz - 1;
  *mc = m;
  *nc = n;
  *pc = p;
}

int extract_data(int operation, double xyz[3], const struct data *in,
                 struct data *out)
{
  int m = 0, n = 0, p = 0;
  closest_index(in, xyz[0], xyz[1], xyz[2], &m, &n, &p);
  if(operation == 1) {
    out->nx = 1;
    out->ny = in->ny;
    out->nz = in->nz;
    out->xmin = xyz[0];
    out->xmax = xyz[0];
    out->ymin = in->ymin;
    out->ymax = in->ymax;
    out->zmin = in->zmin;
    out->zmax = in->zmax;
    int N = in->ny * in->nz;
    out->values = (double*)malloc(N * sizeof(double));
    int i = 0;
    for(p = 0; p < in->nz; p++) {
      for(n = 0; n < in->ny; n++) {
        out->values[i++] = get(in, m, n, p);
      }
    }
  }
  else if(operation == 2) {
    out->nx = in->nx;
    out->ny = 1;
    out->nz = in->nz;
    out->xmin = in->xmin;
    out->xmax = in->xmax;
    out->ymin = xyz[1];
    out->ymax = xyz[1];
    out->zmin = in->zmin;
    out->zmax = in->zmax;
    int N = in->nx * in->nz;
    out->values = (double*)malloc(N * sizeof(double));
    int i = 0;
    for(p = 0; p < in->nz; p++) {
      for(m = 0; m < in->nx; m++) {
        out->values[i++] = get(in, m, n, p);
      }
    }
  }
  else if(operation == 3) {
    out->nx = in->nx;
    out->ny = in->ny;
    out->nz = 1;
    out->xmin = in->xmin;
    out->xmax = in->xmax;
    out->ymin = in->ymin;
    out->ymax = in->ymax;
    out->zmin = xyz[2];
    out->zmax = xyz[2];
    int N = in->nx * in->ny;
    out->values = (double*)malloc(N * sizeof(double));
    int i = 0;
    for(n = 0; n < in->ny; n++) {
      for(m = 0; m < in->nx; m++) {
        out->values[i++] = get(in, m, n, p);
      }
    }
  }
  else if(operation == 4) {
    out->nx = 1;
    out->ny = 1;
    out->nz = 1;
    out->xmin = xyz[0];
    out->xmax = xyz[0];
    out->ymin = xyz[1];
    out->ymax = xyz[1];
    out->zmin = xyz[2];
    out->zmax = xyz[2];
    out->values = (double*)malloc(sizeof(double));
    out->values[0] = get(in, m, n, p);
  }
  else {
    printf("Unknown data extraction operation '%d'\n", operation);
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  if(argc < 2) {
    printf("Usage: %s [option] file(s)\n", argv[0]);
    printf("Option can be\n");
    printf("  -cutx valx: to extract data on plane x = valx\n");
    printf("  -cuty valy: to extract data on plane y = valy\n");
    printf("  -cutz valz: to extract data on plane z = valz\n");
    printf("  -point valx valy valz: to extract data on point (valx, valy, valz)\n");
    printf("If no option is given, extract the full data\n\n");
    printf("The output .msh file can be opened with Gmsh (https://gmsh.info)\n");
    return 0;
  }
  int operation = 0, step = 0;
  double xyz[3] = {0., 0., 0.};
  for(int i = 1; i < argc; i++) {
    if(!strcmp(argv[i], "-cutx") && !step  && i < argc - 1) {
      operation = 1;
      xyz[0] = atof(argv[++i]);
    }
    else if(!strcmp(argv[i], "-cuty") && !step  && i < argc - 1) {
      operation = 2;
      xyz[1] = atof(argv[++i]);
    }
    else if(!strcmp(argv[i], "-cutz") && !step  && i < argc - 1) {
      operation = 3;
      xyz[2] = atof(argv[++i]);
    }
    else if(!strcmp(argv[i], "-point") && !step  && i < argc - 3) {
      operation = 4;
      xyz[0] = atof(argv[++i]);
      xyz[1] = atof(argv[++i]);
      xyz[2] = atof(argv[++i]);
    }
    else {
      struct data d;
      if(read_data(&d, argv[i])) exit(1);
      if(!operation) { // full data
        if(write_data_gmsh(&d, "full.msh", step)) exit(1);
      }
      else {
        struct data e;
        if(extract_data(operation, xyz, &d, &e)) exit(1);
        if(operation == 4) { // point evaluation
          if(write_data_csv(&e, "point.csv", step)) exit(1);
        }
        else { // cut plane
          if(write_data_gmsh(&e, (operation == 1) ? "cutx.msh" :
                             (operation == 2) ? "cuty.msh" : "cutz.msh", step))
            exit(1);
        }
        free_data(&e);
      }
      free_data(&d);
      step++;
    }
  }

  printf("Done!\n");
  return 0;
}
