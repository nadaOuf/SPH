/** File:		sph_system.h
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include "sph_type.h"
#include "slef_def.h"

class Particle
{
public:
	uint id;
	float3 pos;
	float3 vel;
	float3 acc;
	float3 ev;

	float dens;
	float pres;

	float surf_norm;
	////////////////////////***temp****///////////////////////
	Status state;   
	float  temp;           // melting members
	float  temp_eval;
	float heat_fusion;
	float3 particle_color;
	///////////////////////////////////////////////
	Particle *next;
	void CalcParticleColor();
	//Status state;
};

class SPHSystem
{
public:
	uint max_particle;
	uint num_particle;

	float kernel;
	float mass;

	float3 world_size;
	float cell_size;
	uint3 grid_size;
	uint tot_cell;


	float3 IceVelocity;
	float3 IceDeltPos;

	float3 IceForce_fluid;
	float3 IceForce_rigid;
	float3 IceForce_pre;
	float3 IceForce_vis;
	float3 IceForce;


	float3 gravity;
	float wall_damping;
	float rest_density;
	float gas_constant;
	float viscosity;
	float time_step;
	float surf_norm;
	float surf_coe;

	float poly6_value;
	float spiky_value;
	float visco_value;

	float grad_poly6;
	float lplc_poly6;

	float kernel_2;
	float self_dens;
	float self_lplc_color;

	int N_ice;
	int N_w;

	Particle *mem;
	Particle testSource;
	Particle *sPhoton;
	Particle **cell;

	uint sys_running;

	uint Nice;

public:
	SPHSystem();
	~SPHSystem();
	void animation();
	void init_system();
	void add_particle(float3 pos, float3 vel, float T);
	void add_heatSource(float3 pos, float T);

private:
	void build_table();
	void comp_dens_pres();
	void comp_force_adv();
	void advection();
	void compute_Torque();
	void HeatAdvect(Particle *p);
	//void SetColor();
	float HeatTransfer_particle(Particle *pi, Particle *pj);
	//void HeatTransfer();
	float HeatTransferAir(Particle *p, float dA);

private:
	int3 calc_cell_pos(float3 p);
	uint calc_cell_hash(int3 cell_pos);
};

#endif
