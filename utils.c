#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

DATA *new_data(int nx, int ny, int nz, DATA *rho, int initialize){
  DATA *data = malloc(sizeof(DATA));
  if(!data)
  return NULL;

  data->nx = nx;
  data->ny = ny;
  data->nz = nz;
  data->xmin = rho->xmin;
  data->xmax = rho->xmax;
  data->ymin = rho->ymin;
  data->ymax = rho->ymax;
  data->zmin = rho->zmin;
  data->zmax = rho->zmax;

  // values
  data->values = malloc(data->nx*sizeof(double));
  if(!data->values){
    free_data(data);
    return NULL;
  }
  for(int m = 0; m < data->nx; m++){
    data->values[m] = malloc(data->ny*sizeof(double));
    if(!data->values[m]){
      free_data(data);
      return NULL;
    }
    for(int n = 0; n < data->ny; n++){
      data->values[m][n] = malloc(data->nz*sizeof(double));
      if(!data->values[m][n]){
        free_data(data);
        return NULL;
      }
    }
  }
  if (initialize>=0){
    for(int p = 0; p < data->nz; p++){
      for(int n = 0; n < data->ny; n++){
        for(int m = 0; m < data->nx; m++){
          data->values[m][n][p] = initialize;
        }
      }
    }
  }
  return data;
}

DATA *read_data(char *filename){
  DATA *data = malloc(sizeof(DATA));
  if(!data)
  return NULL;
  FILE *fp = fopen(filename,"r");
  if(!fp){
    printf("File couls not be oppened \n");
    free_data(data);
    return NULL;

  }
  int tmp_int;
  double tmp_double;
  // 3 first ints
  if(!fread(&tmp_int,sizeof(int),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->nx = tmp_int;
  if(!fread(&tmp_int,sizeof(int),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->ny = tmp_int;
  if(!fread(&tmp_int,sizeof(int),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->nz = tmp_int;
  // 6 doubles
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->xmin = tmp_double;
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->xmax = tmp_double;
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->ymin = tmp_double;
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->ymax = tmp_double;
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->zmin = tmp_double;
  if(!fread(&tmp_double,sizeof(double),1,fp)){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  data->zmax = tmp_double;
  // values
  data->values = malloc(data->nx*sizeof(double));
  if(!data->values){
    free_data(data);
    fclose(fp);
    return NULL;
  }
  for(int m = 0; m < data->nx; m++){
    data->values[m] = malloc(data->ny*sizeof(double));
    if(!data->values[m]){
      free_data(data);
      fclose(fp);
      return NULL;
    }
    for(int n = 0; n < data->ny; n++){
      data->values[m][n] = malloc(data->nz*sizeof(double));
      if(!data->values[m][n]){
        free_data(data);
        fclose(fp);
        return NULL;
      }
    }
  }
  for(int p = 0; p < data->nz; p++){
    for(int n = 0; n < data->ny; n++){
      for(int m = 0; m < data->nx; m++){
        if(!fread(&tmp_double,sizeof(double),1,fp)){
          free_data(data);
          fclose(fp);
          return NULL;
        }
        data->values[m][n][p] = tmp_double;
      }
    }
  }
  fclose(fp);
  return data;
}

void print_data(DATA *data){
  for(int p = 0; p < data->nz; p++){
    for(int n = 0; n < data->ny; n++){
      for(int m = 0; m < data->nx; m++){
        printf("%f \n",data->values[m][n][p]);
      }
      printf("\n");
    }
  }

}

void write_data(char *filename, DATA *data){
  FILE *fp = fopen(filename,"w");
  if(!fp)
  return;
  // 3 first ints
  fwrite(&data->nx,sizeof(int),1,fp);
  fwrite(&data->ny,sizeof(int),1,fp);
  fwrite(&data->nz,sizeof(int),1,fp);
  // 6 doubles
  fwrite(&data->xmin,sizeof(double),1,fp);
  fwrite(&data->xmax,sizeof(double),1,fp);
  fwrite(&data->ymin,sizeof(double),1,fp);
  fwrite(&data->ymax,sizeof(double),1,fp);
  fwrite(&data->zmin,sizeof(double),1,fp);
  fwrite(&data->zmax,sizeof(double),1,fp);
  // values
  for(int p = 0; p < data->nz; p++){
    for(int n = 0; n < data->ny; n++){
      for(int m = 0; m < data->nx; m++){
        fwrite(&data->values[m][n][p],sizeof(double),1,fp);
      }
    }
  }
  fclose(fp);
  return;
}

void write_data_fromvalues(char *filename, DATA *data, double *values){
  FILE *fp = fopen(filename,"w");
  if(!fp)
  return;
  // 3 first ints
  fwrite(&data->nx,sizeof(int),1,fp);
  fwrite(&data->ny,sizeof(int),1,fp);
  fwrite(&data->nz,sizeof(int),1,fp);
  // 6 doubles
  fwrite(&data->xmin,sizeof(double),1,fp);
  fwrite(&data->xmax,sizeof(double),1,fp);
  fwrite(&data->ymin,sizeof(double),1,fp);
  fwrite(&data->ymax,sizeof(double),1,fp);
  fwrite(&data->zmin,sizeof(double),1,fp);
  fwrite(&data->zmax,sizeof(double),1,fp);
  // values
  int idx =0;
  for(int m = 0; m < data->nx; m++){
    for(int n = 0; n < data->ny; n++){
      for(int p = 0; p < data->nz; p++){
        fwrite(&values[idx],sizeof(double),1,fp);
        idx++;
      }
    }
  }
  fclose(fp);
  return;
}

void swap_data(DATA*** data1, DATA*** data2) {
    DATA** tmp_dat = *data2;
    *data2 = *data1;
    *data1 = tmp_dat;
}


void allocate_data_arrays(DATA*** data1, DATA*** data2, int dimx, int dimy, int dimz, DATA* rho) {
    *data1 = malloc(sizeof(DATA*) * 4);
    *data2 = malloc(sizeof(DATA*) * 4);

    (*data1)[0] = new_data(dimx - 1, dimy, dimz, rho, 0);
    (*data1)[1] = new_data(dimx, dimy - 1, dimz, rho, 0);
    (*data1)[2] = new_data(dimx, dimy, dimz - 1, rho, 0);
    (*data1)[3] = new_data(dimx, dimy, dimz, rho, 0);

    (*data2)[0] = new_data(dimx - 1, dimy, dimz, rho, -1);
    (*data2)[1] = new_data(dimx, dimy - 1, dimz, rho, -1);
    (*data2)[2] = new_data(dimx, dimy, dimz - 1, rho, -1);
    (*data2)[3] = new_data(dimx, dimy, dimz, rho, -1);
}

void free_data(DATA *data){
  for(int i=0;i<data->nx;i++){
    for(int j=0;j<data->ny;j++){
        free(data->values[i][j]);
    }
    free(data->values[i]);
  }
  free(data->values);
  free(data);
  return;
}

void free_data_array(DATA **data, int n){
  for(int i=0;i<n;i++){
    free_data(data[i]);
  }
  free(data);
}

PARAM *read_param(char *filename){
  FILE *fp = fopen(filename,"r");
  if(!fp)
  return NULL;
  PARAM *param = malloc(sizeof(PARAM));
  if(!param){
    fclose(fp);
    return NULL;
  }
  char buffer[500];
  if(!fscanf(fp,"%lf",&param->delta)){
    free(param);
    fclose(fp);
    return NULL;
  }
  if(!fscanf(fp,"%lf",&param->delta_t)){
    free(param);
    fclose(fp);
    return NULL;
  }
  if(!fscanf(fp,"%lf",&param->max_t)){
    free(param);
    fclose(fp);
    return NULL;
  }
  if(!fscanf(fp,"%d",&param->sampling_rate)){
    free(param);
    fclose(fp);
    return NULL;
  }

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->source_type = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->source_type){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->source_type,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->input_speed_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->input_speed_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->input_speed_filename,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->input_density_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->input_density_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->input_density_filename,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->output_pressure_base_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->output_pressure_base_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->output_pressure_base_filename,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->output_velocity_x_base_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->output_velocity_x_base_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->output_velocity_x_base_filename,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->output_velocity_y_base_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->output_velocity_y_base_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->output_velocity_y_base_filename,buffer);

  if(!fscanf(fp,"%s",buffer)){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  param->output_velocity_z_base_filename = malloc((strlen(buffer)+1)*sizeof(char));
  if (!param->output_velocity_z_base_filename){
    free_param(param);
    fclose(fp);
    return NULL;
  }
  strcpy(param->output_velocity_z_base_filename,buffer);

  fclose(fp);
  return param;
}

void free_param(PARAM *param){
  free(param->source_type);
  free(param->input_speed_filename);
  free(param->input_density_filename);
  free(param->output_pressure_base_filename);
  free(param->output_velocity_x_base_filename);
  free(param->output_velocity_y_base_filename);
  free(param->output_velocity_z_base_filename);
  free(param);
  return;
}

BufferData initialize_buffers(int dimx, int dimy, int dimz) {
    BufferData buffers;

    int totalSize = dimx * dimy * dimz;

    buffers.p = (double *) malloc(totalSize * sizeof(double));
    buffers.vx = (double *) malloc(totalSize * sizeof(double));
    buffers.vy = (double *) malloc(totalSize * sizeof(double));
    buffers.vz = (double *) malloc(totalSize * sizeof(double));

    return buffers;
}

void fill_buffer(DATA** data1, DATA** data2, DomainBounds bounds, BufferData *buffers) {
    int buffer_index = 0;

    for (int x = bounds.start_x; x < bounds.end_x; x++) {
      for (int y = bounds.start_y; y < bounds.end_y; y++) {
        for (int z = bounds.start_z; z < bounds.end_z; z++) {
          buffers->p[buffer_index] = data2[3]->values[x][y][z];
          if (x < data1[0]->nx) {
              buffers->vx[buffer_index] = data2[0]->values[x][y][z];
          }
          if (y < data1[1]->ny) {
              buffers->vy[buffer_index] = data2[1]->values[x][y][z];
          }
          if (z < data1[2]->nz) {
              buffers->vz[buffer_index] = data2[2]->values[x][y][z];
          }
          buffer_index++;
        }
      }
    }
}

void free_buffer_data(BufferData *buffers) {
    if (buffers->p) free(buffers->p);
    if (buffers->vx) free(buffers->vx);
    if (buffers->vy) free(buffers->vy);
    if (buffers->vz) free(buffers->vz);
}

// return value of the field x,y,z or 0 if bad index
double get_val(DATA *data,int x,int y,int z){
  if((x < data->nx && x >= 0)&&
     (y < data->ny && y >= 0)&&
     (z < data->nz && z >= 0)){
       return data->values[x][y][z];
     }
  return 0;
}

int map_index(int ni_index, int ni_length, int n_length) {
    return (ni_index*n_length)/ni_length;
}

void updatePressure(DomainBounds bounds, DATA **result,DATA **data, DATA *rho, DATA *speed, PARAM *param){
  double delta_t = param->delta_t;
  double delta = param->delta;
  double diffx;
  double diffy;
  double diffz;
  double coeff;

  // #pragma omp parallel for collapse(3)
      for(int p = bounds.start_z; p < bounds.end_z; p++){
        for(int n = bounds.start_y; n < bounds.end_y; n++){
            for(int m = bounds.start_x; m < bounds.end_x; m++){

        // get_val will return 0 if out of the grid
        diffz = get_val(data[2],m,n,p) - get_val(data[2],m,n,p-1);
        diffy = get_val(data[1],m,n,p) - get_val(data[1],m,n-1,p);
        diffx = get_val(data[0],m,n,p) - get_val(data[0],m-1,n,p);

        int rhom = map_index(m, data[3]->nx, rho->nx);
        int rhon = map_index(n, data[3]->ny, rho->ny);
        int rhop = map_index(p, data[3]->nz, rho->nz);

        int speedm = map_index(m, data[3]->nx, speed->nx);
        int speedn = map_index(n, data[3]->ny, speed->ny);
        int speedp = map_index(p, data[3]->nz, speed->nz);

        coeff = (rho->values[rhom][rhon][rhop] * speed->values[speedm][speedn][speedp] * speed->values[speedm][speedn][speedp] * delta_t / delta);

        result[3]->values[m][n][p] = data[3]->values[m][n][p] - (coeff * (diffx + diffy + diffz));

      }
    }
  }
}

void updateSpeed(DomainBounds bounds, DATA **result,DATA **data, DATA *rho, DATA *speed, PARAM *param){
  double delta_t = param->delta_t;
  double delta = param->delta;
  // Update of speed
  double diffxpressure;
  double diffypressure;
  double diffzpressure;
  double coeff;
  // #pragma omp parallel for collapse(3)
    for(int p = bounds.start_z; p < bounds.end_z; p++){
      for(int n = bounds.start_y; n < bounds.end_y; n++){
        for(int m = bounds.start_x; m < bounds.end_x; m++){

        int rhom = map_index(m, data[3]->nx, rho->nx);
        int rhon = map_index(n, data[3]->ny, rho->ny);
        int rhop = map_index(p, data[3]->nz, rho->nz);

        coeff = (delta_t/(delta*rho->values[rhom][rhon][rhop]));
        
        if(p < data[3]->nz-1){
          diffzpressure = result[3]->values[m][n][p+1]-result[3]->values[m][n][p];
          result[2]->values[m][n][p] = data[2]->values[m][n][p] - (coeff * diffzpressure);
        }
        if(n < data[3]->ny-1){
          diffypressure = result[3]->values[m][n+1][p]-result[3]->values[m][n][p];
          result[1]->values[m][n][p] = data[1]->values[m][n][p] - (coeff * diffypressure);
        }
        if(m < data[3]->nx-1){
          diffxpressure = result[3]->values[m+1][n][p]-result[3]->values[m][n][p];
          result[0]->values[m][n][p] = data[0]->values[m][n][p] - (coeff * diffxpressure);
        }
      }
    }
  }
  for(int i = 0;i<4;i++){
    result[i]->xmin = rho->xmin;
    result[i]->ymin = rho->ymin;
    result[i]->zmin = rho->zmin;
    result[i]->xmax = rho->xmax;
    result[i]->ymax = rho->ymax;
    result[i]->zmax = rho->zmax;
  }
}