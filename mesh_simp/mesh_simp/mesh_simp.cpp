// mesh_simp.cpp : 定义控制台应用程序的入口点。
#include "stdafx.h"
using namespace SimpleOBJ;

static std::vector<Vertex *> vertices;
static std::vector<Triangle *> triangles;
static std::vector<Edge *> edges;

void AddVertex(Vec3f* vert, int n);
void AddFaces(Array<int,3>* tri, int n);

void SimplifyMesh(CSimpleObject &obj, double r);

void ReadFromSimpleObject(const CSimpleObject &obj);
void WriteToSimpleObject(CSimpleObject& obj);

void InitMesh();
void Collapse(Edge * e, int& d); 

//void CompactMesh();

int main(int argc, char* argv[])
{
	assert(argc==4);
	char* input_file=argv[1];
	char* output_file=argv[2];
	double rate=strtod(argv[3], NULL);
	
	/*char *input_file="dinosaur.obj";
	char *output_file="dinosaur_new.obj";
	double rate=0.1;*/
	
	clock_t start=clock();
	
	CSimpleObject obj;
	obj.LoadFromObj(input_file);
 
	SimplifyMesh(obj, rate);

	std::cout<<"Vertex Number = "<<vertices.size()<<std::endl;
	std::cout<<"Triangle Number = "<<triangles.size()<<std::endl;
	obj.SaveToObj(output_file);
	
	clock_t finish=clock();
	double totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	std::cout<<totaltime<<" seconds elapsed."<<std::endl;
	
	getchar();
	return 0;
}

bool cmp(Edge *e1, Edge *e2){
	return e1->cost > e2->cost;
}

void SimplifyMesh(CSimpleObject &obj, double r){
	ReadFromSimpleObject(obj);
	
	InitMesh();
	//for(int i=0 ; i<edges.size(); i++)
	//	std::cout<<edges[i]->cost<<std::endl;
	int face_to_delete = (int) triangles.size() * (1-r);
	int face_deleted = 0;

	while( face_deleted < face_to_delete ){
		pop_heap(edges.begin(), edges.end(), cmp);
		Edge *e= *(edges.end()-1);
		edges.pop_back();
		
		// e may be deleted in previous process
		if (e->deleted)
			continue;
		Collapse(e, face_deleted);
	}
	std::cout<<"Collapse complete."<<std::endl;
	
	for(std::vector<Vertex *>::iterator iter=vertices.begin(); iter!=vertices.end(); ){
		if((*iter)->deleted)
			vertices.erase(iter);
		else
			iter++;
	}
	std::cout<<"Compact complete."<<std::endl;
	
	WriteToSimpleObject(obj);
}

void Collapse(Edge *e, int& d){
	Vertex *u = e->vertex[0];
	Vertex *v = e->vertex[1];

	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			v->Q[i][j]=e->Q[i][j];
	v->position = e->vbar;

	/*for(std::vector<Edge *>::iterator iter=edges.begin();iter!=edges.end();iter++){
		if ((*iter)->deleted)
			continue;
		if((*iter)->vertex[0]==u || (*iter)->vertex[1]==u)
			assert(find(u->pairs.begin(),u->pairs.end(),*iter)!=u->pairs.end();
	}*/

	//delete the two faces.
	for(std::vector<Triangle *>::iterator iter=u->faces.begin(); iter!=u->faces.end(); ){
		if((*iter)->vertex[0]==v||(*iter)->vertex[1]==v||(*iter)->vertex[2]==v){
			triangles.erase(find(triangles.begin(),triangles.end(), *iter));
			delete(*iter);
			d++;
		}
		else
			iter++;
	}

	//refresh neighbors of v.
	for(std::vector<Edge *>::iterator iter=v->pairs.begin(); iter!=v->pairs.end();){
		(*iter)->deleted=true;
		Vertex *u;
		if ((*iter)->vertex[0]==v)
			u=(*iter)->vertex[1];
		else
			u=(*iter)->vertex[0];
		u->pairs.erase(find(u->pairs.begin(), u->pairs.end(), *iter));
		v->pairs.erase(iter);
	}

	for(std::vector<Vertex *>::iterator iter=v->neighbors.begin(); iter!=v->neighbors.end();iter++){
		Edge *e = new Edge(*iter, v);
		//e->ComputeMatrix();
		for(int j=0; j<4; j++)
			for(int k=0; k<4; k++)
				e->Q[j][k]=(*iter)->Q[j][k]+v->Q[j][k];
		e->ComputeVbar();
		e->ComputeCost();
				
		v->pairs.push_back(e);
		(*iter)->pairs.push_back(e);
		edges.push_back(e);
		push_heap(edges.begin(), edges.end(),cmp);
	}

	//for other faces of u, replace u with v.
	for(std::vector<Triangle *>::iterator iter=u->faces.begin(); iter!=u->faces.end(); ){
		if(u == (*iter)->vertex[0])
			(*iter)->vertex[0]=v;
		else if(u == (*iter)->vertex[1])
			(*iter)->vertex[1]=v;
		else{
			assert(u == (*iter)->vertex[2]);
			(*iter)->vertex[2]=v;
		}
		//(*iter)->ComputeNormalMatrix();
		Vertex *tmp[3];
		for(int i=0; i<3; i++)
			tmp[i]=(*iter)->vertex[i];

		v->faces.push_back(*iter);
		Triangle *tmp2= *iter;
		
		u->faces.erase(iter);
		
		//add edge between v and the other two
		for(int i=0;i<3;i++){
			int i2=(i+1)%3;
			if(find(tmp[i]->neighbors.begin(),tmp[i]->neighbors.end(),tmp[i2])==tmp[i]->neighbors.end()){
				tmp[i]->neighbors.push_back(tmp[i2]);
				tmp[i2]->neighbors.push_back(tmp[i]);
				
				Edge *e = new Edge(tmp[i], tmp[i2]);
				//e->ComputeMatrix();
				for(int j=0; j<4; j++)
					for(int k=0; k<4; k++)
						e->Q[j][k]=tmp[i]->Q[j][k]+tmp[i2]->Q[j][k];
				e->ComputeVbar();
				e->ComputeCost();
				
				tmp[i]->pairs.push_back(e);
				tmp[i2]->pairs.push_back(e);
				edges.push_back(e);
				push_heap(edges.begin(), edges.end(),cmp);

			}
		}

		for(int i=0;i<3;i++){
			u->CleanNeighbor(tmp2->vertex[i]);
			tmp2->vertex[i]->CleanNeighbor(u);
		}
	}

	//for now u has no neighbors
	assert(u->neighbors.size()==0);
	assert(u->pairs.size()==0);
	
	/*for(std::vector<Vertex *>::iterator iter=u->neighbors.begin(); iter!=u->neighbors.end();){
		if(find((*iter)->neighbors.begin(), (*iter)->neighbors.end(), u)==(*iter)->neighbors.end())
			std::cout<<"not find "<<d<<std::endl;
		(*iter)->neighbors.erase(find((*iter)->neighbors.begin(), (*iter)->neighbors.end(), u));
		u->neighbors.erase(iter);
	}
		
	for(std::vector<Edge *>::iterator iter=u->pairs.begin(); iter!=u->pairs.end();){
		Vertex * v;
		if ((*iter)->vertex[0]== u)
			v=(*iter)->vertex[1];
		else
			v=(*iter)->vertex[0];
		(*iter)->deleted=true;
		v->pairs.erase(find(v->pairs.begin(), v->pairs.end(), *iter));
		u->pairs.erase(iter);
	}*/
		
	u->deleted=true;
	//delete u;
	//u = NULL;
}

void InitMesh(){
	for(std::vector<Vertex *>::iterator iter=vertices.begin(); iter!=vertices.end(); iter++)
		(*iter)->ComputeMatrix();
	for(std::vector<Edge *>::iterator iter=edges.begin(); iter!=edges.end(); iter++){
		(*iter)->ComputeMatrix();
		(*iter)->ComputeVbar();
		(*iter)->ComputeCost();
	}
	make_heap(edges.begin(), edges.end(), cmp);
	std::cout<<"Construction complete."<<std::endl;
}

void ReadFromSimpleObject(const CSimpleObject &obj){
	AddVertex(obj.m_pVertexList, obj.m_nVertices);
	AddFaces(obj.m_pTriangleList, obj.m_nTriangles);
}

void AddVertex(Vec3f* vert, int n){
	for(int i=0; i<n; i++){
		Vertex *v = new Vertex(vert[i]);
		v->index = i;
		vertices.push_back(v);
	}
	std::cout<<"Vertices Loaded"<<std::endl;
}

void AddFaces(Array<int,3>* tri, int n){
	for(int i=0; i<n; i++){
		Triangle *t=new Triangle(vertices[tri[i][0]],
			vertices[tri[i][1]],
			vertices[tri[i][2]]);
		triangles.push_back(t);
		for (int j=0; j<3; j++){
			t->vertex[j]->faces.push_back(t);
			int k=(j+1)%3;
			if(find(t->vertex[j]->neighbors.begin(), 
				t->vertex[j]->neighbors.end(), 
				t->vertex[k]) == t->vertex[j]->neighbors.end()){
					Edge *e=new Edge(t->vertex[j], t->vertex[k]);
					t->vertex[j]->pairs.push_back(e);
					t->vertex[k]->pairs.push_back(e);
					edges.push_back(e);
					t->vertex[j]->neighbors.push_back(t->vertex[k]);
					t->vertex[k]->neighbors.push_back(t->vertex[j]);
			}
		}
	}
	std::cout<<"Triangles Loaded"<<std::endl;
}

void WriteToSimpleObject(CSimpleObject &obj){
	obj.m_nVertices=vertices.size();
	obj.m_nTriangles=triangles.size();
	
	obj.m_pVertexList = new Vec3f[obj.m_nVertices];
    obj.m_pTriangleList = new Array<int,3> [obj.m_nTriangles];
		
	for(int i=0;i<obj.m_nVertices;i++){
		Vec3f vP(vertices[i]->position);
		obj.m_pVertexList[i] = vP;
	}

    for(int i=0;i<obj.m_nTriangles;i++)
    {
		Array<int,3> vIndices;
		vIndices[0]=find(vertices.begin(), vertices.end(), triangles[i]->vertex[0])-vertices.begin();
		vIndices[1]=find(vertices.begin(), vertices.end(), triangles[i]->vertex[1])-vertices.begin();
		vIndices[2]=find(vertices.begin(), vertices.end(), triangles[i]->vertex[2])-vertices.begin();
		obj.m_pTriangleList[i] = vIndices;
    }
}