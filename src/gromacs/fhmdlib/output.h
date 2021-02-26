#ifndef FHMD_OUTPUT_H_
#define FHMD_OUTPUT_H_

void fhmd_dump_all(FHMD *fh);
void fhmd_write_tecplot_data(FHMD *fh, int step, double time);
void fhmd_AFM_dump_all(FHMD *fh);
void fhmd_shear_dump (FHMD *fh);
void AFM_dump (FHMD *fh);


#endif /* FHMD_OUTPUT_H_ */
