#include"AuxMeshSimp.h"
#include "stdafx.h"
namespace SimpleOBJ
{
	Triangle::Triangle(Vertex *v0, Vertex *v1, Vertex *v2){
		vertex[0]=v0;
		vertex[1]=v1;
		vertex[2]=v2;
		ComputeNormalMatrix();
	}

	Triangle::~Triangle(){
		for(int i=0;i<3;i++){
			//if (vertex[i])
			if (!vertex[i]->deleted)
				vertex[i]->faces.erase(find(vertex[i]->faces.begin(), vertex[i]->faces.end(), this));
		}
		for(int i=0;i<3;i++){
			int i2=(i+1)%3;
			//if(!vertex[i]||!vertex[i2])
			if (vertex[i]->deleted || vertex[i2]->deleted)
				continue;
			vertex[i]->CleanNeighbor(vertex[i2]);
			vertex[i2]->CleanNeighbor(vertex[i]);
		}
	}

	Vec3f Cross(Vec3f u, Vec3f v){
		return Vec3f(u[1]*v[2]-u[2]*v[1], 
			-u[0]*v[2]+u[2]*v[0], 
			u[0]*v[1]-u[1]*v[0]);
	}

	void Triangle::ComputeNormalMatrix(){
		Vec3f v0=vertex[0]->position;
		Vec3f v1=vertex[1]->position;
		Vec3f v2=vertex[2]->position;
		normal = Cross(v1-v0, v2-v0);
		normal.Normalize();

		double temp[4]={ normal.x, normal.y, normal.z,
			-(normal.x * vertex[0]->position.x 
			+ normal.y * vertex[0]->position.y 
			+ normal.z * vertex[0]->position.z)};
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				M[i][j] = temp[i] * temp[j];
			}
		}
	}

	Vertex::Vertex(Vec3f v){
		position = v;
		deleted = false;
		for(int i=0; i<4; i++)
			for (int j=0; j<4; j++)
				Q[i][j]=0;
	}

	void Vertex::ComputeMatrix(){
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				for(std::vector<Triangle *>::iterator iter=faces.begin(); iter!=faces.end(); iter++)
					Q[i][j] += (*iter)-> M[i][j];
	}

	Edge::Edge(Vertex *v0, Vertex *v1){
		vertex[0]=v0;
		vertex[1]=v1;
		deleted=false;
	}

	Edge::~Edge(){
	}

	double Det(double a, double b, double c,
		double d, double e, double f, 
		double g, double h, double i){
			return a*e*i + b*f*g + c*d*h - c*e*g - f*h*a - i*b*d;
	}

	double Edge::AuxComputeCost(Vec3f v){
		double t[4]={v.x, v.y, v.z, 1};
		
		double temp1[4]={0};
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				temp1[i]+=t[j]*Q[j][i];
			}
		}
		double result=0;
		for(int i=0; i<4; i++){
			result+=temp1[i] * t[i];
		}
		return result;
	}
	
	bool Edge::operator == (const Edge other) const{
			return (vertex[0] == other.vertex[0] && vertex[1] == other.vertex[1])||
				(vertex[0] == other.vertex[1] && vertex[1] == other.vertex[0]);
	}

	void Edge::ComputeMatrix(){
		// ax+by+cz+d=0
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				Q[i][j]=vertex[0]->Q[i][j]+vertex[1]->Q[i][j]; 
			}
		}
	}

	void Edge::ComputeVbar(){
		
		double det=Det(Q[0][0], Q[0][1], Q[0][2], 
			Q[1][0], Q[1][1], Q[1][2],
			Q[2][0], Q[2][1], Q[2][2]);

		if(det>1e-9 || det<-1e-9){
			double x= -1/det * Det(Q[0][1], Q[0][2], Q[0][3],
				Q[1][1], Q[1][2], Q[1][3],
				Q[2][1], Q[2][2], Q[2][3]);
			double y= 1/det * Det(Q[0][0], Q[0][2], Q[0][3],
				Q[1][0], Q[1][2], Q[1][3],
				Q[2][0], Q[2][2], Q[2][3]);
			double z= -1/det * Det(Q[0][0], Q[0][1], Q[0][3],
				Q[1][0], Q[1][1], Q[1][3],
				Q[2][0], Q[2][1], Q[2][3]);
			vbar = Vec3f (x, y, z);
			cost = AuxComputeCost(vbar);
			return;
		}
	
	
		double error1=AuxComputeCost(vertex[0]->position);
		double error2=AuxComputeCost(vertex[1]->position);
		double error3=AuxComputeCost((vertex[0]->position+vertex[1]->position)/2);
		cost = error1 < error2 ? (error1 < error3 ? error1: error3):(error2 < error3 ? error2 : error3);
		if(cost == error1) 
			vbar = vertex[0]->position;
		else if(cost == error2)
			vbar = vertex[1]->position;
		else
			vbar= (vertex[0]->position+vertex[1]->position)/2;
		return;
	}
	void Edge::ComputeCost(){
		cost = AuxComputeCost(vbar);
	}

	Vertex::~Vertex(){
		//assert(faces.size()==0);
		for(std::vector<Vertex *>::iterator iter=neighbors.begin(); iter!=neighbors.end();){
			(*iter)->neighbors.erase(find((*iter)->neighbors.begin(), (*iter)->neighbors.end(), this));
			neighbors.erase(iter);
		}
		
		for(std::vector<Edge *>::iterator iter=pairs.begin(); iter!=pairs.end();){
			Vertex * u;
			if ((*iter)->vertex[0]==this)
				u=(*iter)->vertex[1];
			else
				u=(*iter)->vertex[0];
			u->pairs.erase(find(u->pairs.begin(), u->pairs.end(), *iter));
			(*iter)->deleted=true;
			pairs.erase(iter);
		}
	}
	
	void Vertex::CleanNeighbor(Vertex *u){
		if(find(neighbors.begin(),neighbors.end(),u)==neighbors.end())
			return;
		for(int i=0; i<faces.size(); i++)
			if (faces[i]->vertex[0]==u || faces[i]->vertex[1]==u || faces[i]->vertex[2]==u)
				return;
		neighbors.erase(find(neighbors.begin(),neighbors.end(), u));
		for(std::vector<Edge *>::iterator iter=pairs.begin();iter!=pairs.end(); iter++)
			if ((*iter)->vertex[0]==u||(*iter)->vertex[1]==u){
				(*iter)->deleted=true;
				pairs.erase(iter);
				break;
			}
	}
}