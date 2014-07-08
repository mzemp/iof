/*
** iof_tipsy.h
**
** Various handling functions for tipsy format
**
** Written by Marcel Zemp
*/

#ifndef IOF_TIPSY_H
#define IOF_TIPSY_H

#include <rpc/types.h>
#include <rpc/xdr.h>

/*
** Structures
*/

typedef struct tipsy_header {

	double time;
	unsigned int ntotal;
	unsigned int ndim;
	unsigned int ngas;
	unsigned int ndark;
	unsigned int nstar;
	} TIPSY_HEADER;

typedef struct tipsy_gas_particle {

	float mass;
	float pos[3];
	float vel[3];
	float rho;
	float temp;
	float hsmooth;
	float metals;
	float phi;
	} TIPSY_GAS_PARTICLE;

typedef struct tipsy_dark_particle {

	float mass;
	float pos[3];
	float vel[3];
	float eps;
	float phi;
	} TIPSY_DARK_PARTICLE;

typedef struct tipsy_star_particle {

	float mass;
	float pos[3];
	float vel[3];
	float metals;
	float tform;
	float eps;
	float phi;
	} TIPSY_STAR_PARTICLE;

typedef struct tipsy_gas_particle_dpp {

	float mass;
	double pos[3];
	float vel[3];
	float rho;
	float temp;
	float hsmooth;
	float metals;
	float phi;
	} TIPSY_GAS_PARTICLE_DPP;

typedef struct tipsy_dark_particle_dpp {

	float mass;
	double pos[3];
	float vel[3];
	float eps;
	float phi;
	} TIPSY_DARK_PARTICLE_DPP;

typedef struct tipsy_star_particle_dpp {

	float mass;
	double pos[3];
	float vel[3];
	float metals;
	float tform;
	float eps;
	float phi;
	} TIPSY_STAR_PARTICLE_DPP;

typedef struct tipsy_structure {

	TIPSY_HEADER *th;
	TIPSY_GAS_PARTICLE *tgp;
	TIPSY_DARK_PARTICLE *tdp;
	TIPSY_STAR_PARTICLE *tsp;
	} TIPSY_STRUCTURE;

typedef struct tipsy_structure_dpp {

	TIPSY_HEADER *th;
	TIPSY_GAS_PARTICLE_DPP *tgpdpp;
	TIPSY_DARK_PARTICLE_DPP *tdpdpp;
	TIPSY_STAR_PARTICLE_DPP *tspdpp;
	} TIPSY_STRUCTURE_DPP;

/*
** Functions
*/

void copy_tgp_to_tgpdpp(const TIPSY_GAS_PARTICLE *, TIPSY_GAS_PARTICLE_DPP *);
void copy_tdp_to_tdpdpp(const TIPSY_DARK_PARTICLE *, TIPSY_DARK_PARTICLE_DPP *);
void copy_tsp_to_tspdpp(const TIPSY_STAR_PARTICLE *, TIPSY_STAR_PARTICLE_DPP *);
void copy_tgpdpp_to_tgp(const TIPSY_GAS_PARTICLE_DPP *, TIPSY_GAS_PARTICLE *);
void copy_tdpdpp_to_tdp(const TIPSY_DARK_PARTICLE_DPP *, TIPSY_DARK_PARTICLE *);
void copy_tspdpp_to_tsp(const TIPSY_STAR_PARTICLE_DPP *, TIPSY_STAR_PARTICLE *);

void read_tipsy_nb_header(FILE *, TIPSY_HEADER *);
void read_tipsy_nb_gas(FILE *, TIPSY_GAS_PARTICLE *);
void read_tipsy_nb_dark(FILE *, TIPSY_DARK_PARTICLE *);
void read_tipsy_nb_star(FILE *, TIPSY_STAR_PARTICLE *);
void read_tipsy_nb_gas_dpp(FILE *, TIPSY_GAS_PARTICLE_DPP *);
void read_tipsy_nb_dark_dpp(FILE *, TIPSY_DARK_PARTICLE_DPP *);
void read_tipsy_nb_star_dpp(FILE *, TIPSY_STAR_PARTICLE_DPP *);
void read_tipsy_nb(FILE *, TIPSY_STRUCTURE *);
void read_tipsy_nb_dpp(FILE *, TIPSY_STRUCTURE_DPP *);

void write_tipsy_nb_header(FILE *, const TIPSY_HEADER *);
void write_tipsy_nb_gas(FILE *, const TIPSY_GAS_PARTICLE *);
void write_tipsy_nb_dark(FILE *, const TIPSY_DARK_PARTICLE *);
void write_tipsy_nb_star(FILE *, const TIPSY_STAR_PARTICLE *);
void write_tipsy_nb_gas_dpp(FILE *, const TIPSY_GAS_PARTICLE_DPP *);
void write_tipsy_nb_dark_dpp(FILE *, const TIPSY_DARK_PARTICLE_DPP *);
void write_tipsy_nb_star_dpp(FILE *, const TIPSY_STAR_PARTICLE_DPP *);
void write_tipsy_nb(FILE *, const TIPSY_STRUCTURE *);
void write_tipsy_nb_dpp(FILE *, const TIPSY_STRUCTURE_DPP *);

void read_tipsy_xdr_header(XDR *, TIPSY_HEADER *);
void read_tipsy_xdr_gas(XDR *, TIPSY_GAS_PARTICLE *);
void read_tipsy_xdr_dark(XDR *, TIPSY_DARK_PARTICLE *);
void read_tipsy_xdr_star(XDR *, TIPSY_STAR_PARTICLE *);
void read_tipsy_xdr_gas_dpp(XDR *, TIPSY_GAS_PARTICLE_DPP *);
void read_tipsy_xdr_dark_dpp(XDR *, TIPSY_DARK_PARTICLE_DPP *);
void read_tipsy_xdr_star_dpp(XDR *, TIPSY_STAR_PARTICLE_DPP *);
void read_tipsy_xdr(FILE *, TIPSY_STRUCTURE *);
void read_tipsy_xdr_dpp(FILE *, TIPSY_STRUCTURE_DPP *);

void write_tipsy_xdr_header(XDR *, TIPSY_HEADER *);
void write_tipsy_xdr_gas(XDR *, TIPSY_GAS_PARTICLE *);
void write_tipsy_xdr_dark(XDR *, TIPSY_DARK_PARTICLE *);
void write_tipsy_xdr_star(XDR *, TIPSY_STAR_PARTICLE *);
void write_tipsy_xdr_gas_dpp(XDR *, TIPSY_GAS_PARTICLE_DPP *);
void write_tipsy_xdr_dark_dpp(XDR *, TIPSY_DARK_PARTICLE_DPP *);
void write_tipsy_xdr_star_dpp(XDR *, TIPSY_STAR_PARTICLE_DPP *);
void write_tipsy_xdr(FILE *, const TIPSY_STRUCTURE *);
void write_tipsy_xdr_dpp(FILE *, const TIPSY_STRUCTURE_DPP *);

void read_tipsy_ascii(FILE *, TIPSY_STRUCTURE *);
void read_tipsy_ascii_dpp(FILE *, TIPSY_STRUCTURE_DPP *);

void write_tipsy_ascii(FILE *, const TIPSY_STRUCTURE *);
void write_tipsy_ascii_dpp(FILE *, const TIPSY_STRUCTURE_DPP *);

#endif /* IOF_TIPSY_H */
