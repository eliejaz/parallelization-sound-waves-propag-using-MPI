#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

void InttoString(unsigned long integer,unsigned long length,char *out){
    unsigned long rest = 0;
    unsigned long div = (unsigned long)pow(10,length-2);

    for(int x=0;x<length-1;x++){
      out[x] = (integer-rest)/div;
      rest = rest + (out[x]*div);
      div = div/10;
      out[x]=out[x]+48;
    }
    out[length-1] = '\0';
    return;
}

int get_domain_size_x(int rank, int size, int dimx) {
    int subdomain_size_x = dimx / size;
    int remainder = dimx % size;
    int start_x, end_x;

    if (rank < remainder) {
        // This process gets an extra element
        start_x = rank * (subdomain_size_x + 1);
        end_x = start_x + subdomain_size_x + 1;
    } else {
        // This process gets the regular number of elements
        start_x = rank * subdomain_size_x + remainder; 
        end_x = start_x + subdomain_size_x;
    }

    return end_x - start_x;
}


void write_data_fromvalues_mpi(char *filename, DATA *data, double *values, int rank, int size){
    MPI_File fh;
    MPI_Status status;
    MPI_Offset offset;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Rank 0 will write the headers
    if (rank == 0) {
        MPI_File_write(fh, &data->nx, 1, MPI_INT, &status);
        MPI_File_write(fh, &data->ny, 1, MPI_INT, &status);
        MPI_File_write(fh, &data->nz, 1, MPI_INT, &status);
        MPI_File_write(fh, &data->xmin, 1, MPI_DOUBLE, &status);
        MPI_File_write(fh, &data->xmax, 1, MPI_DOUBLE, &status);
        MPI_File_write(fh, &data->ymin, 1, MPI_DOUBLE, &status);
        MPI_File_write(fh, &data->ymax, 1, MPI_DOUBLE, &status);
        MPI_File_write(fh, &data->zmin, 1, MPI_DOUBLE, &status);
        MPI_File_write(fh, &data->zmax, 1, MPI_DOUBLE, &status);
    }
    // Calculate the offset for each rank in bytes
    int header_size = 3 * sizeof(int) + 6 * sizeof(double);
    int data_size_per_rank_bytes = get_domain_size_x(rank, size, data->nx) * data->ny * data->nz * sizeof(double);
    
    offset = header_size;
    for(int i = 0; i < rank; i++) {
        offset += get_domain_size_x(i, size, data->nx) * data->ny * data->nz * sizeof(double);
    }
    MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, values, data_size_per_rank_bytes, MPI_BYTE, &status);
    MPI_File_close(&fh);
}

void write_simulation_data(int iteration, int int_size, PARAM* p, DATA** data2,  BufferData *buffers, int rank, int size) {
    char str_int[int_size];
    char buffer[50];

    InttoString(iteration, int_size, str_int);

    // Write velocity x data
    strcpy(buffer, p->output_velocity_x_base_filename);
    write_data_fromvalues_mpi(strcat(buffer, str_int), data2[0], buffers->vx, rank, size);

    // Write velocity y data
    strcpy(buffer, p->output_velocity_y_base_filename);
    write_data_fromvalues_mpi(strcat(buffer, str_int), data2[1], buffers->vy, rank, size);

    // Write velocity z data
    strcpy(buffer, p->output_velocity_z_base_filename);
    write_data_fromvalues_mpi(strcat(buffer, str_int), data2[2], buffers->vz, rank, size);

    // Write pressure data
    strcpy(buffer, p->output_pressure_base_filename);
    write_data_fromvalues_mpi(strcat(buffer, str_int), data2[3], buffers->p, rank, size);
}

void check_arguments(int argc) {
    if (argc != 2) {
        printf("Error with parameters\nExpected: ./a.out file_name.txt\n");
        exit(-1);
    }
}

void read_simulation_data(char* param_name, PARAM** p, DATA** rho, DATA** speed) {
    *p = read_param(param_name);
    *rho = read_data((*p)->input_density_filename);
    *speed = read_data((*p)->input_speed_filename);
}

void calculate_dimensions(PARAM* p, DATA* rho, int* dimx, int* dimy, int* dimz) {
    *dimx = (rho->xmax - rho->xmin) / p->delta;
    *dimy = (rho->ymax - rho->ymin) / p->delta;
    *dimz = (rho->zmax - rho->zmin) / p->delta;

    if(*dimz == 0){
        *dimz = 1;
    }
}

DomainBounds calculate_domain_bounds(int rank, int size, int dimx, int dimy, int dimz) {
    DomainBounds bounds;

    int subdomain_size_x = dimx / size;
    int remainder = dimx % size;
    if (rank < remainder) {
        // This process gets an extra element
        bounds.start_x = rank * (subdomain_size_x + 1);
        bounds.end_x = bounds.start_x + subdomain_size_x + 1;
    } else {
        // This process gets the regular number of elements
        bounds.start_x = rank * subdomain_size_x + remainder; 
        bounds.end_x = bounds.start_x + subdomain_size_x;
    }
    printf("rank: %i -- startx: %i -- endx: %i \n", rank, bounds.start_x, bounds.end_x);


    bounds.start_y = 0;
    bounds.end_y = dimy;
    printf("rank: %i -- starty: %i -- endy: %i \n", rank, bounds.start_y, bounds.end_y);

    bounds.start_z = 0;
    bounds.end_z = dimz;
    printf("rank: %i -- startz: %i -- endz: %i \n", rank, bounds.start_z, bounds.end_z);


    return bounds;
}

void send_receive_data_p(int rank, int size, DomainBounds bounds, DATA** data2, double* sendReceiveBuffer) {
    int dimy = bounds.end_y - bounds.start_y;
    int dimz = bounds.end_z - bounds.start_z;
    int bufSize = dimy * dimz;
    MPI_Request send_request, recv_request;

    if (rank > 0) {
        for (int y = bounds.start_y; y < bounds.end_y; y++) {
            for (int z = bounds.start_z; z < bounds.end_z; z++) {
                sendReceiveBuffer[y * dimz + z] = data2[3]->values[bounds.start_x][y][z];
            }
        }
        MPI_Isend(sendReceiveBuffer, bufSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_request);
    }
    if (rank < size - 1) {
        MPI_Irecv(sendReceiveBuffer, bufSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv_request);
    }

    // Wait for non-blocking operations to complete
    if (rank > 0) {
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
        for (int y = bounds.start_y; y < bounds.end_y; y++) {
            for (int z = bounds.start_z; z < bounds.end_z; z++) {
                data2[3]->values[bounds.end_x][y][z] = sendReceiveBuffer[y * dimz + z];
            }
        }
    }
}


void send_receive_data_v(int rank, int size, DomainBounds bounds, DATA** data2, double* sendReceiveBuffer) {
    int dimy = bounds.end_y - bounds.start_y;
    int dimz = bounds.end_z - bounds.start_z;
    int bufSize = dimy * dimz;
    MPI_Request send_request, recv_request;

    if (rank > 0) {
        MPI_Irecv(sendReceiveBuffer, bufSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv_request);
    }
    if (rank < size - 1) {
        for (int y = bounds.start_y; y < bounds.end_y; y++) {
            for (int z = bounds.start_z; z < bounds.end_z; z++) {
                sendReceiveBuffer[y * dimz + z] = data2[0]->values[bounds.end_x - 1][y][z];
            }
        }
        MPI_Isend(sendReceiveBuffer, bufSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_request);
    }

    // Wait for non-blocking operations to complete
    if (rank > 0) {
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
        for (int y = bounds.start_y; y < bounds.end_y; y++) {
            for (int z = bounds.start_z; z < bounds.end_z; z++) {
                data2[0]->values[bounds.start_x - 1][y][z] = sendReceiveBuffer[y * dimz + z];
            }
        }
    }
    if (rank < size - 1) {
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
    }
}

void fill_original_data(DATA** data1, int dimx, int dimy, int dimz, int i, PARAM* p) {
    if (strcmp(p->source_type, "point_source_middle_3400") == 0) {
        data1[3]->values[dimx / 2][dimy / 2][0] = sin(2 * M_PI * 3400 * i * p->delta_t);
    } 
    else if (strcmp(p->source_type, "point_source_middle_6800") == 0) {
        data1[3]->values[dimx / 2][dimy / 2][0] = sin(2 * M_PI * 6800 * i * p->delta_t);
    } 
    else if (strcmp(p->source_type, "uniform_source") == 0) {
        data1[3]->values[dimx / 2][dimy / 2][0] = 1.0;
    }
    else if (strcmp(p->source_type, "oscillating_source_1000") == 0) {
        data1[3]->values[dimx / 2][dimy / 2][0] = sin(2 * M_PI * 1000 * i * p->delta_t);
    }
    else if (strcmp(p->source_type, "random_noise_source") == 0) {
        data1[3]->values[dimx / 2][dimy / 2][0] = (double)rand() / RAND_MAX;
    }
    else {
        // ToDo
    }
}

void performSimulationStep(int rank, int size, DomainBounds bounds, DATA** data1, DATA** data2,
                           DATA* rho, DATA* speed, PARAM* p, double* sendReceivedBuffer) {
    updatePressure(bounds, data2, data1, rho, speed, p);
    send_receive_data_p(rank, size, bounds, data2, sendReceivedBuffer);
    updateSpeed(bounds, data2, data1, rho, speed, p);
    send_receive_data_v(rank, size, bounds, data2, sendReceivedBuffer);
}

void run_simulation(int rank, int size, PARAM* p, DATA** data1, DATA** data2, DATA* rho, DATA* speed, int dimx, int dimy, int dimz) {

    int bufSize = dimy * dimz;
    double *sendRecievedbuffer = (double *) malloc(bufSize * sizeof(double));

    DomainBounds bounds = calculate_domain_bounds(rank, size, dimx, dimy, dimz);
    BufferData buffers = initialize_buffers(dimx, dimy, dimz);

    for (int i = 0; i < (int)(p->max_t / p->delta_t); i++) {
      fill_original_data(data1, dimx, dimy, dimz, i, p);
      performSimulationStep(rank, size, bounds, data1, data2, rho, speed, p, sendRecievedbuffer);
      
      if (!((i + 1) % p->sampling_rate)) {
         fill_buffer(data1, data2, bounds, buffers);
        if(rank==0){
          write_simulation_data(i, 9, p, data2, &buffers, rank, size);
        }
      }
      swap_data(&data1, &data2);
      MPI_Barrier(MPI_COMM_WORLD);
    }

    free(sendRecievedbuffer);
    free_buffer_data(&buffers);

}

int main(int argc, char *argv[]) {
    check_arguments(argc);

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PARAM *p;
    DATA *rho, *speed;
    read_simulation_data(argv[1], &p, &rho, &speed);

    int dimx, dimy, dimz;
    calculate_dimensions(p, rho, &dimx, &dimy, &dimz);

    DATA **data1, **data2;
    allocate_data_arrays(&data1, &data2, dimx, dimy, dimz, rho);

    run_simulation(rank, size, p, data1, data2, rho, speed, dimx, dimy, dimz);

    free_data_array(data1, 4);
    free_data_array(data2, 4);
    free_param(p);
    free_data(rho);
    free_data(speed);

    MPI_Finalize();
    return 0;
}