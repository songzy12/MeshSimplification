#pragma once
#include "SimpleObject.h"
#include <vector>
namespace SimpleOBJ
{
	class Vertex;
	class Edge;
	class Triangle{
	public:
		Vertex *vertex[3];
		
		Vec3f normal;
		double M[4][4];
		
		Triangle(Vertex *v0, Vertex *v1, Vertex *v2);
		~Triangle();
		
		void ComputeNormalMatrix();
	};
	
	class Vertex{
	public:
		int index;
		Vec3f position;
		bool deleted;

		std::vector<Edge *> pairs;
		std::vector<Triangle *> faces;
		std::vector<Vertex *> neighbors;

		double Q[4][4];
		void ComputeMatrix();
		void CleanNeighbor(Vertex *u);

		Vertex(Vec3f v);
		~Vertex();
	};

	class Edge{
	public:
		Vertex *vertex[2];
		
		double Q[4][4];
		double cost;
		Vec3f vbar;
		bool deleted;

		double AuxComputeCost(Vec3f v);
		void ComputeMatrix();
		void ComputeCost();
		void ComputeVbar();
		bool operator == (const Edge other) const;
		Edge(Vertex *v0, Vertex *v1);
		~Edge();
	};
}