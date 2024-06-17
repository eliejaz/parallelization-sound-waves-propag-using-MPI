#ifndef UTILS_H
#define UTILS_H

typedef struct {
  int nx,ny,nz; //number of nodes along x,y,z ; they must all be > 0
  double xmin,xmax,ymin,ymax,zmin,zmax; // min/max x,y,z coordinate of the domain wtf is that ??
  double ***values; //values along the grid (nx,ny,nz)
} DATA;

typedef struct {
  double delta; // spatial grid offset
  double delta_t; // time step
  double max_t; // max simulation time
  int sampling_rate;
  char *source_type;
  char *input_speed_filename;
  char *input_density_filename;
  char *output_pressure_base_filename;
  char *output_velocity_x_base_filename;
  char *output_velocity_y_base_filename;
  char *output_velocity_z_base_filename;
} PARAM;

typedef struct {
    int start_x;
    int end_x;
    int start_y;
    int end_y;
    int start_z;
    int end_z;
} DomainBounds;

typedef struct {
    double *p;
    double *vx;
    double *vy;
    double *vz;
} BufferData;

DATA *new_data(int nx, int ny, int nz, DATA *rho, int initialize);

DATA *read_data(char *filename);

void print_data(DATA *data);

void write_data(char *filename, DATA *data);

void write_data_fromvalues(char *filename, DATA *data, double *values);

void swap_data(DATA*** data1, DATA*** data2);

void allocate_data_arrays(DATA*** data1, DATA*** data2, int dimx, int dimy, int dimz, DATA* rho);

void free_data(DATA *data);

void free_data_array(DATA **data, int n);

PARAM *read_param(char *filename);

void free_param(PARAM *param);

BufferData initialize_buffers(int dimx, int dimy, int dimz);

void fill_buffer(DATA** data1, DATA** data2, DomainBounds bounds, BufferData *buffers);

void free_buffer_data(BufferData *buffers);

double get_val(DATA *data,int x,int y,int z);

int map_index(int ni_index, int ni_length, int n_length);

void updatePressure(DomainBounds bounds, DATA **result,DATA **data, DATA *rho, DATA *speed, PARAM *param);

void updateSpeed (DomainBounds bounds, DATA **result,DATA **data, DATA *rho, DATA *speed, PARAM *param);

void update_step(DATA **result,DATA **data, DATA *rho, DATA *speed, PARAM *param);

#endif
