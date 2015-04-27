/** File:		sph_system.cpp
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

#include "sph_system.h"
#include "sph_header.h"

#include "Partio.h"
#include <sstream>
#include<iostream>
using namespace std;
SPHSystem::SPHSystem()
{
	frame_number = 0;
	max_particle=30000;
	num_particle=0;

	kernel=0.04f;
	mass=0.02f;
	//mass=0.0008f;
	world_size.x=0.64f;
	world_size.y=0.64f;
	world_size.z=0.64f;
	cell_size= kernel;
	grid_size.x=(uint)ceil(world_size.x/cell_size);
	grid_size.y=(uint)ceil(world_size.y/cell_size);
	grid_size.z=(uint)ceil(world_size.z/cell_size);
	tot_cell=grid_size.x*grid_size.y*grid_size.z;

	gravity.x=0.0f; 
	gravity.y=-9.8f;
	gravity.z=0.0f;
	wall_damping=-0.5f;
	rest_density=1000.0f;
	gas_constant=1.0f;
	viscosity=3.7f;
	time_step= 0.003f;
	surf_norm=6.0f;
	surf_coe=0.1f;

	poly6_value=315.0f/(64.0f * PI * pow(kernel, 9));;
	spiky_value=-45.0f/(PI * pow(kernel, 6));
	visco_value=45.0f/(PI * pow(kernel, 6));

	grad_poly6=-945/(32 * PI * pow(kernel, 9));
	lplc_poly6=-945/(8 * PI * pow(kernel, 9));

	kernel_2=kernel*kernel;
	self_dens=mass*poly6_value*pow(kernel, 6);
	self_lplc_color=lplc_poly6*mass*kernel_2*(0-3/4*kernel_2);

	mem=(Particle *)malloc(sizeof(Particle)*max_particle);
	cell=(Particle **)malloc(sizeof(Particle *)*tot_cell);

	sys_running=0;

	printf("Initialize SPH:\n");
	printf("World Width : %f\n", world_size.x);
	printf("World Height: %f\n", world_size.y);
	printf("World Length: %f\n", world_size.z);
	printf("Cell Size  : %f\n", cell_size);
	printf("Grid Width : %u\n", grid_size.x);
	printf("Grid Height: %u\n", grid_size.y);
	printf("Grid Length: %u\n", grid_size.z);
	printf("Total Cell : %u\n", tot_cell);
	printf("Poly6 Kernel: %f\n", poly6_value);
	printf("Spiky Kernel: %f\n", spiky_value);
	printf("Visco Kernel: %f\n", visco_value);
	printf("Self Density: %f\n", self_dens);
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
}

int tt = 0;

void SPHSystem::animation()
{

	//tt++;
	if(sys_running == 0)
	{
		return;
	}
	++frame_number;

	float3 pos,vel;
	//pos.x = 0;
	pos.y = world_size.y;
	//pos.z = 0;
	vel.x = 0;
	vel.y = -1;
	vel.z = 0;
	/*
	if(tt%3==0)
	for(pos.x=world_size.x*0.00f; pos.x<world_size.x*0.15f; pos.x+=(kernel*0.7f))
	{
		//for(pos.y=world_size.y*0.00f; pos.y<world_size.y*0.1f; pos.y+=(kernel*0.5f))
		{
			for(pos.z=world_size.z*0.00f; pos.z<world_size.z*0.15f; pos.z+=(kernel*0.7f))
			{
				add_particle(pos, vel , 800);
			}
		}
	}*/

	build_table();
	comp_dens_pres();
	comp_force_adv();
	advection();

}

void SPHSystem::init_system(std::string file_path)
{
	float3 pos;
	float3 vel;

	vel.x=0.0f;
	vel.y=0.0f;
	vel.z=0.0f;

	pos.x=0;//world_size.x;
	pos.y=world_size.y/2;
	pos.z=0;//world_size.z;

	add_heatSource(pos, 400);

	//for(pos.x=world_size.x*0.2f; pos.x<world_size.x*0.3f; pos.x+=(kernel*0.5f))
	//	for(pos.y=world_size.y*0.3f; pos.y<world_size.y*0.8f; pos.y+=(kernel*0.5f))
	//	{
	//		for(pos.z=world_size.z*0.2f; pos.z<world_size.z*0.4f; pos.z+=(kernel*0.5f))
	//		{
	//			add_particle(pos, vel, 200);
	//		}
	//	}
	//for(pos.x=world_size.x*0.3f+kernel*0.5f; pos.x<world_size.x*0.5f; pos.x+=(kernel*0.5f))//x*0.3,x*0.4
	//		for(pos.z=world_size.z*0.2f; pos.z<world_size.z*0.4f; pos.z+=(kernel*0.5f))
	//		{
	//			for(pos.y=world_size.y*0.5f; pos.y<world_size.y*0.6f; pos.y+=(kernel*0.5f))
	//	{
	//			add_particle(pos, vel,200);
	//		}
	//	}
	//for(pos.x=world_size.x*0.5f+kernel*0.5f; pos.x<world_size.x*0.65f; pos.x+=(kernel*0.5f))
	//	for(pos.y=world_size.y*0.3f; pos.y<world_size.y*0.8f; pos.y+=(kernel*0.5f))
	//	{
	//		for(pos.z=world_size.z*0.2f; pos.z<world_size.z*0.4f; pos.z+=(kernel*0.5f))
	//		{
	//			add_particle(pos, vel,200);
	//		}
	//	}

	for(pos.x=world_size.x*0.2f; pos.x<world_size.x*0.6f; pos.x+=(kernel*0.5f))
	{
		for(pos.y=world_size.y*0.3f; pos.y<world_size.y*0.9f; pos.y+=(kernel*0.5f))
		{
			for(pos.z=world_size.z*0.2f; pos.z<world_size.z*0.6f; pos.z+=(kernel*0.5f))
			{
				add_particle(pos, vel, 200);
			}
		}
	}

	printf("Init Particle: %u\n", num_particle);

	//Set the output path for the output files
	output_file_path = file_path;
}

void SPHSystem::add_particle(float3 pos, float3 vel, float T)
{
	Particle *p=&(mem[num_particle]);

	p->id=num_particle;

	p->pos=pos;
	p->vel=vel;

	p->acc.x=0.0f;
	p->acc.y=0.0f;
	p->acc.z=0.0f;
	p->ev.x=0.0f;
	p->ev.y=0.0f;
	p->ev.z=0.0f;

	p->dens=rest_density;
	p->pres=0.0f;

	p->next=NULL;

	p->particle_color.x = 255;
	p->particle_color.y = 0;
	p->particle_color.z = 0;

	p->temp = T ;

	if(p->temp<=273)
	{
		N_ice++;
		p->state = SOLID;
	}
	else
	{
		N_w++;
		p->state = LIQUID;
	}

	p->heat_fusion = 0;

	num_particle++;
}

void SPHSystem::add_heatSource(float3 pos, float T)
{
	testSource.pos = pos;
	testSource.temp = T;
}

void SPHSystem::build_table()
{
	Particle *p;
	uint hash;

	for(uint i=0; i<tot_cell; i++)
	{
		cell[i]=NULL;
	}

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);
		hash=calc_cell_hash(calc_cell_pos(p->pos));

		if(cell[hash] == NULL)
		{
			p->next=NULL;
			cell[hash]=p;
		}
		else
		{
			p->next=cell[hash];
			cell[hash]=p;
		}
	}
}

void SPHSystem::comp_dens_pres()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->dens=0.0f;
		p->pres=0.0f;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					while(np != NULL)
					{
						rel_pos.x=np->pos.x-p->pos.x;
						rel_pos.y=np->pos.y-p->pos.y;
						rel_pos.z=np->pos.z-p->pos.z;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2<INF || r2>=kernel_2)
						{
							np=np->next;
							continue;
						}

						p->dens=p->dens + mass * poly6_value * pow(kernel_2-r2, 3);

						np=np->next;
					}
				}
			}
		}

		p->dens=p->dens+self_dens;
		p->pres=(pow(p->dens / rest_density, 7) - 1) *gas_constant;
	}
}

void SPHSystem::comp_force_adv()
{
	Nice = 0;
	IceForce_fluid.x = 0;
	IceForce_fluid.y = 0;
	IceForce_fluid.z = 0;

	IceForce_rigid.x = 0;
	IceForce_rigid.y = 0;
	IceForce_rigid.z = 0;

	IceVelocity.x = 0;
	IceVelocity.y = 0;
	IceVelocity.z = 0;

	IceDeltPos.x = 0;
	IceDeltPos.y = 0;
	IceDeltPos.z = 0;

	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float3 rel_vel;

	float r2;
	float r;
	float kernel_r;
	float V;

	float pres_kernel;
	float visc_kernel;
	float temp_force;

	float3 grad_color;
	float lplc_color;

	bool solid = false;

	Status pState;
	float ni = 0.0f;
	for(uint i=0; i<num_particle; i++)
	{
		ni = 0.0f;
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		pState = p->state;

		p->acc.x=0.0f;
		p->acc.y=0.0f;
		p->acc.z=0.0f;


		////////////////***solid-Boundary Checking***/////////////////////////
		//////////////////////////////////////////////////////////
		if(pState==SOLID)
		{
			Nice++;
			/////////////////bot///////////
			if(p->pos.y < 0.0f)
			{
				IceForce_rigid.y = -gravity.y - p->vel.y*Bounce/time_step;
				IceVelocity.y = p->vel.y*wall_damping;
				IceDeltPos.y = 0.0 - p->pos.y;
			}

			////////////////////////top///////////////
			if(p->pos.y >world_size.y-BOUNDARY){
			    IceForce_rigid.y = -gravity.y - p->vel.y*Bounce/time_step;
				IceVelocity.y = p->vel.y*wall_damping;
				IceDeltPos.y = world_size.y - p->pos.y;
			}
			/////////////////////////////////////right//////////
			if(p->pos.x >world_size.x-BOUNDARY){
			    IceForce_rigid.x = - p->vel.x*Bounce/time_step;
				IceVelocity.x = p->vel.x*wall_damping;
			    IceDeltPos.x = world_size.x - p->pos.x;
			}
			////////////////left//////////
			if(p->pos.x <0.0){
			    IceForce_rigid.x = - p->vel.x*Bounce/time_step;
				IceVelocity.x = p->vel.x*wall_damping;
			    IceDeltPos.x = 0.0 - p->pos.x;
			}
			//////////////////////////front///////
			if(p->pos.z <0.0){
			    IceForce_rigid.z = - p->vel.z*Bounce/time_step;
				IceVelocity.z = p->vel.z*wall_damping;
			    IceDeltPos.z = 0.0 - p->pos.z;
			}

			//////////////back///////////
			if(p->pos.z >world_size.z-BOUNDARY){
			    IceForce_rigid.z = - p->vel.z*Bounce/time_step;
				IceVelocity.z = p->vel.z*wall_damping;
			    IceDeltPos.z = world_size.z - p->pos.z;
			}
			solid = true;
		}
		////////////////////////////////////////////////////////
		//////////////////////////////////
		grad_color.x=0.0f;
		grad_color.y=0.0f;
		grad_color.z=0.0f;
		lplc_color=0.0f;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						//no neighbor particles :
						//1. check boundary/rigid
						continue;
					}					
					np=cell[hash];
					if(np != NULL && (x + y + z != 0) && ((x == 0 && y == 0) || (x == 0 && z == 0) || ( y == 0 && z == 0)))
						++ni; //sum the number of particles surrounding

					if(np==NULL)
					{//Surface - check heat source
						float3 dist;
						dist.x =  testSource.pos.x - p->pos.x;
						dist.y =  testSource.pos.y - p->pos.y;
						dist.z =  testSource.pos.z - p->pos.z;
						float d2 = dist.x*dist.x+dist.y*dist.y+dist.z*dist.z;
						//Calc Source Heat
						//change this to ray trace......slow....???
						if(dist.x*x>=0 && dist.y*y>=0 && dist.z*z>=0)//vector cosin stuff...
						{
							if(p->temp >= 275 && p->state == SOLID)
								p->heat_fusion += EPSILON*(testSource.temp-p->temp)/d2;
							else
							p->temp += EPSILON*(testSource.temp-p->temp)*time_step/(d2);
						}

					}
					while(np!=NULL)
					{
						p->temp_eval+=HeatTransfer_particle(p, np);
						//cout << p->temp_eval << endl;
						//if(pState)
						rel_pos.x=p->pos.x-np->pos.x;
						rel_pos.y=p->pos.y-np->pos.y;
						rel_pos.z=p->pos.z-np->pos.z;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;
						//p->temp_eval+=HeatTransfer_particle(p, np);
						if(r2 < kernel_2 && r2 > INF)
						{
							r=sqrt(r2);
							V=mass/np->dens;
							kernel_r=kernel-r;

							pres_kernel=spiky_value * kernel_r * kernel_r;
							temp_force=V * (p->pres+np->pres) * pres_kernel;

							p->acc.x=p->acc.x-rel_pos.x*temp_force/r;
							p->acc.y=p->acc.y-rel_pos.y*temp_force/r;
							p->acc.z=p->acc.z-rel_pos.z*temp_force/r;
							float pre = temp_force/(r);

							rel_vel.x=np->ev.x-p->ev.x;
							rel_vel.y=np->ev.y-p->ev.y;
							rel_vel.z=np->ev.z-p->ev.z;

							
							visc_kernel=visco_value*(kernel-r);
							temp_force=V * viscosity * visc_kernel;
							p->acc.x=p->acc.x + rel_vel.x*temp_force; 
							p->acc.y=p->acc.y + rel_vel.y*temp_force; 
							p->acc.z=p->acc.z + rel_vel.z*temp_force; 

							float vis = temp_force;
							//vis = 0;
							// Interfacial intension
							if(p->state==SOLID && np->state==LIQUID) 
							{
								IceForce_fluid.x -= 1.4*rel_vel.x*(vis+pre)/np->dens;
								IceForce_fluid.y -= 1*rel_vel.y*(vis+pre)/np->dens;
								IceForce_fluid.z -= 1.4*rel_vel.z*(vis+pre)/np->dens;
							}

							//if((x + y + z != 0) && ((x == 0 && y == 0) || (x == 0 && z == 0) || ( y == 0 && z == 0)))
							//interfacial tension fi.
							if(np->state == LIQUID)

							{
								p->acc.x=p->acc.x - rel_pos.x*(Kw/r2)/mass; 
								p->acc.y=p->acc.y - rel_pos.y*(Kw/r2)/mass; 
								p->acc.z=p->acc.z - rel_pos.z*(Kw/r2)/mass; 
							}
							else if(np->state == SOLID)
							{
								p->acc.x=p->acc.x - rel_pos.x*(Kice/r2)/mass;
								p->acc.y=p->acc.y - rel_pos.y*(Kice/r2)/mass; 
								p->acc.z=p->acc.z - rel_pos.z*(Kice/r2)/mass; 
							}
							

							float temp=(-1) * grad_poly6 * V * pow(kernel_2-r2, 2);
							grad_color.x += temp * rel_pos.x;
							grad_color.y += temp * rel_pos.y;
							grad_color.z += temp * rel_pos.z;
							lplc_color += lplc_poly6 * V * (kernel_2-r2) * (r2-3/4*(kernel_2-r2));
						}
						np=np->next;
					}
					////////
					//p->temp+=p->temp_eval*time_step;
					/////////////
				}
			}
		}

		float dA = (6 - ni) / 6;
		if(ni > 6) 
			dA = 0.0f;

		float cd = 0.0f;
		if(p->state==LIQUID)
			cd = THERMAL_CONDUCTIVITY_WATER;
		if(p->state==SOLID)
			cd = THERMAL_CONDUCTIVITY_ICE;
		if(p->temp >= 275 && p->state == SOLID)
			p->heat_fusion += cd*p->temp_eval;
		else
			p->temp+= cd*p->temp_eval*time_step;

		//Add the heat transfer due to air 
			
		//Determining the thermal conductivity of the particle depending on its state
		if(p->state==LIQUID)
			cd = HEAT_CAPACITY_WATER;
		if(p->state==SOLID)
			cd = HEAT_CAPACITY_ICE;
		if(p->temp >= 275 && p->state == SOLID)
			p->heat_fusion += (HeatTransferAir(p, dA)/(cd*mass));//+ p->temp_eval;
		else
			p->temp += (HeatTransferAir(p, dA)/(cd*mass))*time_step;
		

		lplc_color+=self_lplc_color/p->dens;
		p->surf_norm=sqrt(grad_color.x*grad_color.x+grad_color.y*grad_color.y+grad_color.z*grad_color.z);

		if(p->surf_norm > surf_norm)
		{
			p->acc.x+=surf_coe * lplc_color * grad_color.x / p->surf_norm;
			p->acc.y+=surf_coe * lplc_color * grad_color.y / p->surf_norm;
			p->acc.z+=surf_coe * lplc_color * grad_color.z / p->surf_norm;
		}
	}

	if(!solid)
		viscosity = 0.5f;

}

void SPHSystem::advection()
{
	int temp_Nice = Nice;
	//cout<<temp_Nice<<endl;
	N_ice = 0;
	N_w = 0;
	Particle *p;
	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

		if(p->state == SOLID)
		{
			p->acc.x = (IceForce_fluid.x)/(temp_Nice+1) + (IceForce_rigid.x)*p->dens;
			p->acc.y = (IceForce_fluid.y)/(temp_Nice+1) + (IceForce_rigid.y)*p->dens;
			p->acc.z = (IceForce_fluid.z)/(temp_Nice+1) + (IceForce_rigid.z)*p->dens;
		}
		//latent heat
		if(p->state == SOLID) 
		{
			if(p->temp>=275 && p->heat_fusion > 300)
			{
				p->state = LIQUID;
			}

		}
		if(p->state == LIQUID)
		{
			N_w++;
		}

		p->vel.x=p->vel.x+p->acc.x*time_step/p->dens+gravity.x*time_step;
		p->vel.y=p->vel.y+p->acc.y*time_step/p->dens+gravity.y*time_step;
		p->vel.z=p->vel.z+p->acc.z*time_step/p->dens+gravity.z*time_step;

		p->pos.x=p->pos.x+p->vel.x*time_step;
		p->pos.y=p->pos.y+p->vel.y*time_step;
		p->pos.z=p->pos.z+p->vel.z*time_step;

		if(p->state == SOLID)
		{
			p->pos.x += IceDeltPos.x; 
			p->pos.y += IceDeltPos.y; 
			p->pos.z += IceDeltPos.z; 
			continue;
		}

		if(p->pos.x >= world_size.x-BOUNDARY)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=world_size.x-BOUNDARY;
		}

		if(p->pos.x < 0.0f)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=0.0f;
		}

		if(p->pos.y >= world_size.y-BOUNDARY)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=world_size.y-BOUNDARY;
		}

		if(p->pos.y < 0.0f)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=0.0f;
		}

		if(p->pos.z >= world_size.z-BOUNDARY)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=world_size.z-BOUNDARY;
		}

		if(p->pos.z < 0.0f)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=0.0f;
		}

		p->ev.x=(p->ev.x+p->vel.x)/2;
		p->ev.y=(p->ev.y+p->vel.y)/2;
		p->ev.z=(p->ev.z+p->vel.z)/2;
	}
	p->CalcParticleColor();
}
void SPHSystem::compute_Torque(){
	//IceForce= IceForce_rigid+IceForce_fluid;
//
}
float SPHSystem::HeatTransferAir(Particle *p, float dA)
{
	//First calculate Q
	float Q = 0.0f, cd = 0.0f;
	Q = THERMAL_CONDUCTIVITY_AIR * (T_air - p->temp) * dA;

	//return the temprature change due to air
	return Q;
}
float SPHSystem::HeatTransfer_particle(Particle *pj, Particle *pi){
	float distx,disty,distz;
    float rij;
	float cd;
	 float smooth_k;
	float temp_neighborEffect = 0;
    distx = pj->pos.x - pi->pos.x;
	disty = pj->pos.y - pi->pos.y;
	distz = pj->pos.z - pi->pos.z;

	rij  = sqrt(pow(distx,2)+ pow(disty,2)+ pow(distz,2));

	if(rij == 0)
		return 0;
	if(rij<R_HEATAFFECT){

			smooth_k=(45.0/(PI*pow(R_HEATAFFECT,6)))*(R_HEATAFFECT-rij);
			temp_neighborEffect= -1*mass*(pj->temp-pi->temp)*smooth_k/pj->dens;
			//cout<<temp_neighborEffect << endl;
			//pi->temp_eval+=temp_neighborEffect;
  			return temp_neighborEffect;

         }
	//	else return 0.0;

}
	

//void SPHSystem::_SetColor(){
	/*Particle *p;
	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);
		if(p->temp<ICE_T){
	       p->particle_color.x=0.0;
	       p->particle_color.y=0.0;
	       p->particle_color.z=p->temp/(ICE_T-MIN_T);	
	      }
	    if(p->temp>ICE_T){
	       p->particle_color.x=0.0;
	       p->particle_color.y=0.0;
	       p->particle_color.z=p->temp/(MAX_T-ICE_T);
        }
	}*/
//}
int3 SPHSystem::calc_cell_pos(float3 p)
{
	int3 cell_pos;//0.0043 f
	cell_pos.x = int(floor((p.x) / cell_size));//0.0399
	cell_pos.y = int(floor((p.y) / cell_size));
	cell_pos.z = int(floor((p.z) / cell_size));

    return cell_pos;
}

uint SPHSystem::calc_cell_hash(int3 cell_pos)
{
	if(cell_pos.x<0 || cell_pos.x>=(int)grid_size.x || cell_pos.y<0 || cell_pos.y>=(int)grid_size.y || cell_pos.z<0 || cell_pos.z>=(int)grid_size.z)
	{
		return (uint)0xffffffff;
	}

	cell_pos.x = cell_pos.x & (grid_size.x-1);  
    cell_pos.y = cell_pos.y & (grid_size.y-1);  
	cell_pos.z = cell_pos.z & (grid_size.z-1);  

	return ((uint)(cell_pos.z))*grid_size.y*grid_size.x + ((uint)(cell_pos.y))*grid_size.x + (uint)(cell_pos.x);
}

//Function that generates a file with all particle information
void SPHSystem::partio_file_output()
{
	//If animation is not running yet
	if(sys_running == 0)
	{
		return;
	}

	//build the file name which is frame + frame_number + ext
	std::string frame_str = "frame";

	//Convert from int to string using stringstream
	std::stringstream frame_number_str;
	frame_number_str << frame_number;

	std::string ext_str = ".bgeo";

	std::string file_string = output_file_path + "\\" + frame_str + frame_number_str.str() + ext_str;

	//Create partio object
	Partio::ParticlesDataMutable *partio_data = Partio::create();
	
	//Set the number of particles
	partio_data->addParticles(num_particle);

	//Create attribute varibles
	Partio::ParticleAttribute position_attribute = partio_data->addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute velocity_attribute = partio_data->addAttribute("v", Partio::VECTOR, 3);
	Partio::ParticleAttribute id_attribute		 = partio_data->addAttribute("id", Partio::INT, 1);

	Particle *p;
	//Loop through the particles and set the attributes
	for(int i = 0; i < num_particle; ++i)
	{
		p=&(mem[i]);

		//Set the position in the partio object
		float* pos = partio_data->dataWrite<float>(position_attribute, i);
		pos[0] = p->pos.x *10;
		pos[1] = p->pos.y *10;
		pos[2] = p->pos.z *10;

		//Set the velocity in the partio object
		float* v = partio_data->dataWrite<float>(velocity_attribute, i);
		v[0] = p->vel.x;
		v[1] = p->vel.y;
		v[2] = p->vel.z;

		//Set the id of the particle
		int* id = partio_data->dataWrite<int>(id_attribute, i);
		id[0] = i;
	}

	//Write out the file
	Partio::write(file_string.c_str(), *partio_data);
	
	partio_data->release();	
	
}

void Particle::CalcParticleColor()
{
	float dv = MAX_T - MIN_T;
	float3 RGB;
	RGB.x = 1;
	RGB.y = 1;
	RGB.z = 1;

	if (temp < MIN_T) temp = MIN_T;
	if (temp > MAX_T) temp = MAX_T;
	
	if (temp < ( MIN_T+ 0.25 * dv)) {
		RGB.x = 0;
		RGB.y = 4 * (temp - MIN_T) / dv;
	} 
	else if (temp < (MIN_T + 0.5 * dv)) 
	{
		RGB.x = 0.0;
        RGB.z = 1.0 + 4.0 * (MIN_T + 0.25 * dv - temp) / dv;
	} 
	else if (temp < (MIN_T + 0.75 * dv))
	{
		RGB.x = 4.0 * (temp - MIN_T - 0.5 * dv) / dv;
        RGB.z = 0.0;
	} 
	else {
		RGB.y = 1.0 + 4.0 * (MIN_T + 0.75 * dv - temp) / dv;
		RGB.z = 0.0;
	}
	RGB.x*=255;
	RGB.y*=255;
	RGB.z*=255;
	particle_color.x = RGB.x;
	particle_color.y = RGB.y;
	particle_color.z = RGB.z;
}