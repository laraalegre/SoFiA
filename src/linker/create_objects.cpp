#ifndef Connected_Gegion_H
#define Connected_Gegion_H
#include <vector>
#include <algorithm>
#include "types.h"
using namespace std;

class Similarity
{
public:
	Similarity()
	{id=0;sameas=0;};
	Similarity(int _id)
	{id=_id;sameas=_id;};
	int id, sameas, tag;
};


class ConnectedComponents
{
public:
	ConnectedComponents();
	~ConnectedComponents();
	vector<Similarity> labels;
	int highest_label;
	vector<int>result;
	back_strct bk_st;
	//---------------------------------------------------------------------
	ConnectedComponents(int soft_maxlabels) ;
	void clear();
	bool is_root_label(int id);
	
	int root_of(int id);
	bool is_equivalent(int id, int as);
	bool merge(int id1, int id2);
	int new_label();
	void label_image(double *img,bool K8_connectivity);
	double*imgpt;
	int *labelimg;
	int *label(double *pixp){return &labelimg[pixp-imgpt];}
	double thresh;
	int left,right,top,bottom;
	int above(double v);
	int above(double v,int x,int y);
	int index(double*pixp){return int(pixp-imgpt);}
	int relabel_image();
	int obnum;
	void alloc_space();
	int width;
	int height;
	vector<ObjInfo> ResultVec;
	void GetResultGroup();
	void assign_back_strt();
	bool if_dif_sigma;
};

void ConnectedComponents::assign_back_strt()
{
	bk_st.width = width;
	bk_st.height = height;
	bk_st.nbackx = (bk_st.width -1)/bk_st.backw + 1;
	bk_st.nbacky = (bk_st.height-1)/bk_st.backh + 1;
	bk_st.nback = bk_st.nbackx * bk_st.nbacky;
	bk_st.fvec.resize(bk_st.nback);
}

int ConnectedComponents::above(double v,int x,int y)
{
	int index_x,index_y;
	index_x=x/bk_st.backw;
	index_y=y/bk_st.backh;
	int sig_index;
	sig_index=index_y*bk_st.nbackx+index_x;
	return v > bk_st.fvec[sig_index] ? 1:0;
}

void ConnectedComponents::GetResultGroup()
{
	//----------------------------------------------------
	int i;
	ResultVec.clear();
	ResultVec.resize(obnum-1);
	for(i=0;i<result.size();i++)
	{
		fSPOINT pt;
		pt.x=result[i]%width;
		pt.y=result[i]/width;
		pt.Flux=imgpt[result[i]] ;
		pt.dthresh=thresh;
		int group_num= labelimg[result[i]];
		if(group_num>(obnum-1))
			continue;
		ResultVec[group_num-1].PtVec.push_back(pt);
	}
	//----------------------------------------------------
	vector<ObjInfo> ResultVecBk;
	ResultVecBk.clear();
	for(i=0;i<ResultVec.size();i++)
	{
		ObjInfo Tmp;
		Tmp=ResultVec[i];
		if(Tmp.PtVec.size()>3)
			ResultVecBk.push_back(Tmp);
	}
	ResultVec=ResultVecBk;
}

ConnectedComponents::ConnectedComponents()
{
	labelimg=NULL;
}

ConnectedComponents::~ConnectedComponents()
{
	delete []labelimg;
	labelimg=NULL;
}

void ConnectedComponents::alloc_space()
{
	delete []labelimg;
	labelimg=NULL;
	int cnt;cnt=width*height;
	labelimg=new int[cnt];
	int i;
	for(i=0;i<cnt;i++)labelimg[i]=0;
}

int ConnectedComponents::above(double v)
{
	if(if_dif_sigma)
		return above(v,bk_st.crtx,bk_st.crty);
	return v>thresh? 1:0;
};

bool ConnectedComponents::is_root_label(int id)
{return (labels[id].sameas == id);};

void ConnectedComponents::clear()
{ 
	result.clear();
	fill(labels.begin(), labels.end(), Similarity());
	highest_label = 0;
};

ConnectedComponents::ConnectedComponents(int soft_maxlabels)
: labels(soft_maxlabels)
{clear();};

int ConnectedComponents::root_of(int id)
{
	while (!is_root_label(id)) {
		labels[id].sameas = labels[labels[id].sameas].sameas;
		id = labels[id].sameas;
	}
	return id;
} ;

bool ConnectedComponents::is_equivalent(int id, int as)
{
	return (root_of(id) == root_of(as));
};

bool ConnectedComponents::merge(int id1, int id2)
{
	if(!is_equivalent(id1, id2))
	{
		labels[root_of(id1)].sameas = root_of(id2);
		return false;
	}
	return true;
};

int ConnectedComponents::new_label()
{
	if(unsigned int(highest_label+1) > labels.capacity())
		labels.reserve(highest_label*2);
	labels.resize(highest_label+1);
	labels[highest_label] = Similarity(highest_label);
	return highest_label++;
}


void ConnectedComponents::label_image(double *img,bool K8_CONNECTIVITY)
{
	imgpt=img;
	double *row=0;
	double *last_row = 0;
	int crtval;
	int tempv1,tempv2;
	
	clear();
	row = &img[width*top];
	bk_st.crtx=left;
	bk_st.crty=top;
	if(above(row[left]))
	{
		new_label();
		new_label();
		*label(&row[left])=1;
		result.push_back(index(&row[left]));
	}
	else
	{
		new_label();
		*label(&row[left])=0;
	}
	for(int c=left+1; c<=right; ++c)
	{
		bk_st.crtx=c;
		bk_st.crty=top;
		int crtval=above(row[c]);
		
		if(crtval==0)
		{
			*label(&row[c])=0;
			continue;
		}
		else
		{
			bk_st.crtx=c-1;
			bk_st.crty=top;
			if(above(row[c-1]))
				*label(&row[c]) = *label(&row[c-1]);
			else
				*label(&row[c]) = new_label();
			
			result.push_back(index(&row[c]));
		}
	}
	
	for(int r=top+1; r<=bottom; ++r)
	{
		last_row = row;
		row = &img[width*r];
		
		bk_st.crtx=left;
		bk_st.crty=r;
		crtval=above(row[left]);
		
		if(crtval)
		{
			bk_st.crtx=left;
			bk_st.crty=r-1;
			if(crtval== above(last_row[left]))
				*label(&row[left]) = *label(&last_row[left]);
			else
				*label(&row[left]) = new_label();
			result.push_back(index(&row[left]));
		}
		else
		{
			*label(&row[left])=0;
		}
		for(int c=left+1; c<=right; ++c)
		{
			bk_st.crtx = c;
			bk_st.crty = r;
			crtval=above(row[c]);
			
			bk_st.crtx = c-1;
			bk_st.crty = r;
			tempv1 = above(row[c-1]);
			bk_st.crtx = c;
			bk_st.crty = r-1;
			tempv2=above(last_row[c]);
			
			if(K8_CONNECTIVITY && tempv1 && tempv2)
				merge(*label(&row[c-1]), *label(&last_row[c]));
			
			if(crtval==0)
			{
				*label(&row[c])=0;
				continue;
			}
			
			result.push_back(index(&row[c]));
			
			int mylab = -1;
			
			bk_st.crtx=c-1;
			bk_st.crty=r;
			if(above(row[c-1]))
				mylab = *label(&row[c-1]);
			
			for(int d=(K8_CONNECTIVITY?-1:0); d<1; ++d)
			{
				bk_st.crtx=c+d;
				bk_st.crty=r-1;
				if(above(last_row[c+d]))
				{
					if(mylab>0)
						merge(mylab, *label(&last_row[c+d]));
					else
						mylab = *label(&last_row[c+d]);
				}
			}
			
			if(mylab>0)
				*label(&row[c]) = mylab;
			else
				*label(&row[c]) = new_label();
			
		}
	}
}

int ConnectedComponents::relabel_image()
{
	int newtag = 0;
	int i;
	for(int id=0; id<labels.size(); ++id)
	{if(is_root_label(id))
		labels[id].tag = newtag++;
	}
	
	for(i = 0; i<width*height; i++)
	{
		if(labelimg[i]!=0)
		{
			
			labelimg[i] = labels[root_of(labelimg[i])].tag;
		}
	}
	obnum=newtag;
	return newtag;
}

#endif 