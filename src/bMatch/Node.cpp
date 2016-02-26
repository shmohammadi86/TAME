#include "bMatching.h"
#include "mtxReader.h"
#include "pcl_stack.h"
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <float.h>
#include <string.h>

using namespace std;

int pivotID;
float pivotW;
int ptype;

void Node::print()
{
    int i;
    printf("Printing the heap array\n");
    for(i=0;i<curSize;i++)
    {
        printf("%d %f\n",heap[i].id,heap[i].weight);
    }
}

bool heapComp(Info left, Info right)
{
    return (left.weight<right.weight || (left.weight==right.weight && left.id<right.id));
}


void Node::AddHeap(float wt, int idx )
{
    if(curSize==maxSize)
    {
        if(maxSize>2)
        {
            //heap[0].weight=wt;
            //heap[0].id=idx;

            /// Only heapify one branch of the heap tree
            int small,ri,li,pi=0;
            int done=0;

            if(heap[2].weight >heap[1].weight || (heap[2].weight==heap[1].weight && heap[2].id>heap[1].id))
                small=1;
            else
                small=2;

            if(wt>heap[small].weight || (wt==heap[small].weight && idx>heap[small].id))
            {
                heap[0].weight=heap[small].weight;
                heap[0].id=heap[small].id;
                heap[small].weight=wt;
                heap[small].id=idx;
            }
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }

            pi=small;
            while(!done)
            {
                li=2*pi+1;
                ri=2*pi+2;
                small=pi;

                if(li <maxSize && (heap[li].weight< heap[small].weight || (heap[li].weight == heap[small].weight && heap[li].id < heap[small].id )))
                    small=li;
                if(ri <maxSize && (heap[ri].weight< heap[small].weight || (heap[ri].weight == heap[small].weight && heap[ri].id < heap[small].id)))
                    small=ri;

                if(small != pi)
                {
                    wt=heap[pi].weight;
                    idx=heap[pi].id;

                    heap[pi].weight=heap[small].weight;
                    heap[pi].id=heap[small].id;

                    heap[small].weight=wt;
                    heap[small].id=idx;
                }
                else done=1;

                pi=small;
            }
        }
        else
        {
            if(maxSize != 1 && (wt>heap[1].weight || (wt==heap[1].weight && idx > heap[1].id)))
            {
                heap[0].weight=heap[1].weight;
                heap[0].id=heap[1].id;
                heap[1].weight=wt;
                heap[1].id=idx;
            }
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }
        }

        minEntry.id=heap[0].id;
        minEntry.weight=heap[0].weight;
    }
    else
    {
        heap[curSize].weight=wt;
        heap[curSize].id=idx;
        curSize++;
        if(curSize==maxSize)
        {
            sort(heap,heap+curSize,heapComp);
            minEntry.weight=heap[0].weight;
            minEntry.id=heap[0].id;
        }
    }
}

void Node::Add(float wt, int idx )
{
    if(curSize==maxSize)
    {
        if(maxSize>2)
        {

            //heap[0].weight=wt;
            //heap[0].id=idx;

            int small,ri,li,pi=0;
            int done=0;

            if(heap[2].weight >heap[1].weight || (heap[2].weight==heap[1].weight && heap[2].id>heap[1].id))
                small=1;
            else
                small=2;

            if(wt>heap[small].weight || (wt==heap[small].weight && idx>heap[small].id))
            {
                heap[0].weight=heap[small].weight;
                heap[0].id=heap[small].id;
                heap[small].weight=wt;
                heap[small].id=idx;
            }
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }

            pi=small;
            while(!done)
            {
                li=2*pi+1;
                ri=2*pi+2;
                small=pi;

                if(li <maxSize && (heap[li].weight< heap[small].weight || (heap[li].weight == heap[small].weight && heap[li].id < heap[small].id )))
                    small=li;
                if(ri <maxSize && (heap[ri].weight< heap[small].weight || (heap[ri].weight == heap[small].weight && heap[ri].id < heap[small].id)))
                    small=ri;

                if(small != pi)
                {
                    wt=heap[pi].weight;
                    idx=heap[pi].id;

                    heap[pi].weight=heap[small].weight;
                    heap[pi].id=heap[small].id;

                    heap[small].weight=wt;
                    heap[small].id=idx;
                }
                else done=1;

                pi=small;
            }
        }
        else
        {
            if(maxSize != 1 &&(wt>heap[1].weight || (wt==heap[1].weight && idx > heap[1].id)))
            {
                heap[0].weight=heap[1].weight;
                heap[0].id=heap[1].id;
                heap[1].weight=wt;
                heap[1].id=idx;
            }
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }

        }
        minEntry.id=heap[0].id;
        minEntry.weight=heap[0].weight;
    }
    else
    {
        switch(curSize)
        {
            case 0:
                heap[0].weight=wt;
                heap[0].id=idx;
                break;
            case 1:
                if(heap[0].weight< wt || (heap[0].weight == wt && heap[0].id < idx))
                {
                    heap[1].weight=wt;
                    heap[1].id=idx;
                }
                else
                {
                    heap[1].weight=heap[0].weight;
                    heap[1].id=heap[0].id;

                    heap[0].weight=wt;
                    heap[0].id=idx;
                }
                break;
            case 2:
                if(heap[0].weight< wt || (heap[0].weight == wt && heap[0].id < idx))
                {
                    heap[2].weight=wt;
                    heap[2].id=idx;
                }
                else
                {
                    heap[2].weight=heap[0].weight;
                    heap[2].id=heap[0].id;

                    heap[0].weight=wt;
                    heap[0].id=idx;
                }
                break;
            default:
                heap[curSize].weight=wt;
                heap[curSize].id=idx;

                int pi,ci=curSize;
                int done=0;

                while(!done)
                {
                    pi=(int)floor((ci-1)/2.0);

                    if(heap[pi].weight< heap[ci].weight || (heap[pi].weight == heap[ci].weight && heap[pi].id < heap[ci].id ))
                        done=1;
                    else
                    {
                        wt=heap[pi].weight;
                        idx=heap[pi].id;

                        heap[pi].weight=heap[ci].weight;
                        heap[pi].id=heap[ci].id;

                        heap[ci].weight=wt;
                        heap[ci].id=idx;
                    }

                    if(pi==0)
                        break;

                    ci=pi;
                }
                break;

        }// switch
        curSize++;
        if(curSize==maxSize)
        {
            minEntry.weight=heap[0].weight;
            minEntry.id=heap[0].id;
        }
    }

}

bool comparatorE(EdgeE left, EdgeE right)
{
    return (left.weight > right.weight || (left.weight==right.weight && left.id > right.id));
}

void multiQselectIter(Edge* verInd, int start, int cstart, int cend, int* order, int ostart, int oend, stack_pcl<Param>* st)
{
    int i,p,r,id,length;
    double weight;
    Param ptemp;

    ptemp.cstart=cstart;
    ptemp.cend=cend;
    ptemp.ostart=ostart;
    ptemp.oend=oend;

    while(true)
    {
        if(ptemp.cstart>ptemp.cend || ptemp.ostart > ptemp.oend)
        {
            if(st->empty())
                break;

            ptemp=st->pop();

            //cstart=ptemp.cstart;
            //cend=ptemp.cend;
            //ostart=ptemp.ostart;
            //oend=ptemp.oend;
        }

        p=ptemp.cstart;
        r=ptemp.cend;
        weight = verInd[(p+r)/2].weight;
        id=verInd[(p+r)/2].id;
        while ( p < r )
        {
            while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
                p++;
            while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
                r--;

            if(verInd[p].id==verInd[r].id)
                p++;
            else
                if( p < r )
                {
                    /*temp.id= verInd[p].id;
                    temp.weight=verInd[p].weight;

                    verInd[p].id  = verInd[r].id;
                    verInd[p].weight = verInd[r].weight;

                    verInd[r].id = temp.id;
                    verInd[r].weight = temp.weight;*/

                    swap(verInd[r],verInd[p]);
                }
        }

        length=r-start+1;

        for(i=ptemp.ostart;i<=ptemp.oend;i++)
            if(length<=order[i])
                break;
       cstart=ptemp.cstart;
       ostart=ptemp.ostart;
        if(i<=ptemp.oend)
        {
            if(order[i]==length)
            {
                //multiQselect(verInd,start,r+1,cend,order,i+1,oend);
                ptemp.cstart=r+1;
                //ptemp.cend=cend;
                ptemp.ostart=i+1;
                //ptemp.oend=oend;
                st->push(ptemp);
            }
            else
            {
                //multiQselect(verInd,start,r+1,cend,order,i,oend);

                ptemp.cstart=r+1;
                //ptemp.cend=cend;
                ptemp.ostart=i;
                //ptemp.oend=oend;
                st->push(ptemp);
            }
        }

        //multiQselect(verInd,start,cstart,r-1,order,ostart,i-1);
        ptemp.cstart=cstart;
        ptemp.ostart=ostart;
        ptemp.cend=r-1;
        ptemp.oend=i-1;
    }
}

void multiQselect(Edge* verInd, int start, int cstart, int cend, int* order, int ostart, int oend)
{
    if(cstart>cend || ostart > oend)
        return;

    int i,p,r,id,length;
    double weight;
    Edge temp;

    p=cstart;
    r=cend;
    weight = verInd[(p+r)/2].weight;
    id=verInd[(p+r)/2].id;
    while ( p < r )
    {
        while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
            p++;
        while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
            r--;

        if(verInd[p].id==verInd[r].id)
            p++;
        else
            if( p < r )
            {
                temp.id= verInd[p].id;
                temp.weight=verInd[p].weight;

                verInd[p].id  = verInd[r].id;
                verInd[p].weight = verInd[r].weight;

                verInd[r].id = temp.id;
                verInd[r].weight = temp.weight;
            }
    }

    length=r-start+1;

    for(i=ostart;i<=oend;i++)
        if(length<=order[i])
            break;

    multiQselect(verInd,start,cstart,r-1,order,ostart,i-1);

    if(i<=oend)
    {
        if(order[i]==length)
            multiQselect(verInd,start,r+1,cend,order,i+1,oend);
        else
            multiQselect(verInd,start,r+1,cend,order,i,oend);
    }
}

void my_small_sort(Edge * verInd, int num)
{

#if 0
	for(int size = num; size > 1; size--)
	{
		int size_aligned = 0;
		double min = DBL_MAX; int location = -1;
#ifdef __MIC__
                size_aligned = size/8 * 8;
		for(int i = 0; i < size_aligned; i+=8)
		{
			__m512d v_i = _mm512_undefined_pd();
			v_i = _mm512_loadunpacklo_pd(v_i, (void *)(verInd + i));
			v_i = _mm512_loadunpackhi_pd(v_i, (void *)(verInd + i + 8));

			double local_min = _mm512_reduce_min_pd(v_i);
			if(min > local_min)
			{
				min = local_min;
				__mmask8 v_mask = _mm512_cmpeq_pd_mask(v_i, _mm512_set1_pd(local_min));
				location = i + tzcnt_32(v_mask);
			}
		}
#endif
		for(int i = size_aligned; i < size; i++)
		{
			double v_i = *(double *)(verInd + i);
			if(min > v_i) {
				min = v_i;
				location = i;
			}
		}
		double end = *(double *)(verInd + size-1); *(double *)(verInd+ size-1) = min; *(double *)(verInd + location) = end;
	}
#else
	// insertion sort ... from Wikipedia :)
	for(int i = 1; i < num; i++)
	{
		int j = i;
		while( (j > 0) && ( (verInd[j - 1].weight < verInd[j].weight) || (verInd[j - 1].weight== verInd[j].weight && verInd[j-1].id < verInd[j].id)))
		{
			Edge tmp = verInd[j]; verInd[j] = verInd[j-1]; verInd[j-1] = tmp;
			j--;
		}
	}

#endif
}


int custom_sort_optimized(Edge* verInd, int start, int end, int step, int* order, int type, Edge* neighbors)
{
    int part=start+step;
    int tstart, tend,k,length;
    int id,p,r;
    double weight;
    Edge temp;
    Edge median;

    switch(type)
    {
        case 1: break;
        case 2: part=end;
                sort(verInd+start,verInd+end,comparator);
                //csort(verInd,start,end);
                break;

        case 3: if(part>=end)
                {
                    part=end;
                    sort(verInd+start,verInd+end,comparator);
                }
                else
                {
                    int last=(end-start)/step;
                    //if(last>50)
                        //last=50;
                    //multiQselectIter(verInd,start,start,end-1,order,0,last-1,st);
                    multiQselect(verInd,start,start,end-1,order,0,last-1);
                    sort(verInd+start,verInd+part,comparator);
                }
                break;
        case 4: if(part>=end)
                {
                    part=end;
                    sort(verInd+start,verInd+end,comparator);
                    //my_small_sort(verInd+start, (end  - start));
                }
                else
                {
                    tend=end-1;
                    tstart=start;
                    k=step+1;

                    while(true)
                    {
                        p=tstart;
                        r=tend;
                        weight = verInd[(r+p)/2].weight;
                        id=verInd[(r+p)/2].id;
			median = verInd[(r+p)/2];

#if defined _ACTIVATE_MIC_CODE_ && defined __MIC__
			int num = (tend - tstart + 1);
			int num_aligned = num/16 * 16;
			SIMDFPTYPE v_weight = _MM_SET(weight);
			SIMDINTTYPE v_id = _MM_SET_INT(id);
			Edge * write_1 = verInd + tstart;
			Edge * write_2 = neighbors;
			int cnt_write_1 = 0;
			int cnt_write_2 = 0;

			for(int i = tstart; i < tstart + num_aligned; i+=16)
			{
				// AOS, hence two loads required to get 16 elements
				SIMDINTTYPE v_i = _MM_LOADU_INT((int *)(verInd + i));
				SIMDINTTYPE v_i_16 = _MM_LOADU_INT((int *)(verInd + i + 8));

				// First element in structure is id, and second is weight.
				// Mask with 0xAAAA to get weight, and 0x5555 to get id
				SIMDMASKTYPE v_gt_weight_1 = _MM_CMP_LT_MASK(0xAAAA, v_weight, _mm512_castsi512_ps(v_i));
				SIMDMASKTYPE v_eq_weight_1 = _MM_CMP_EQ_MASK(0xAAAA, _mm512_castsi512_ps(v_i), v_weight);
				SIMDMASKTYPE v_gt_id_1     = _MM_CMP_LT_MASK_INT(0x5555, v_id, v_i);
				SIMDMASKTYPE v_mask_left_1 = _mm512_kor(v_gt_weight_1, _mm512_kand(v_eq_weight_1, (v_gt_id_1<<1)));
				unsigned int cnt_1 = countbits_32(v_mask_left_1);
				v_mask_left_1 = ((v_mask_left_1 >> 1)) | v_mask_left_1;
				_MM_STOREU_MASK_INT((int *)(write_1), v_mask_left_1, v_i);
				write_1 += cnt_1; cnt_write_1 += cnt_1;

				SIMDMASKTYPE v_lt_weight_1 = _MM_CMP_LT_MASK(0xAAAA, _mm512_castsi512_ps(v_i), v_weight);
				SIMDMASKTYPE v_lt_id_1     = _MM_CMP_LT_MASK_INT(0x5555, v_i, v_id);
				SIMDMASKTYPE v_mask_right_1 = _mm512_kor(v_lt_weight_1, _mm512_kand(v_eq_weight_1, (v_lt_id_1<<1)));
				unsigned int cnt_right_1 = countbits_32(v_mask_right_1);
				v_mask_right_1 = ((v_mask_right_1 >> 1)) | v_mask_right_1;
				_MM_STOREU_MASK_INT((int *)(write_2), v_mask_right_1, v_i);
				write_2 += cnt_right_1; cnt_write_2 += cnt_right_1;


				SIMDMASKTYPE v_gt_weight_2 = _MM_CMP_LT_MASK(0xAAAA, v_weight, _mm512_castsi512_ps(v_i_16));
				SIMDMASKTYPE v_eq_weight_2 = _MM_CMP_EQ_MASK(0xAAAA, _mm512_castsi512_ps(v_i_16), v_weight);
				SIMDMASKTYPE v_gt_id_2     = _MM_CMP_LT_MASK_INT(0x5555, v_id, v_i_16);
				SIMDMASKTYPE v_mask_left_2 = _mm512_kor(v_gt_weight_2, _mm512_kand(v_eq_weight_2, (v_gt_id_2<<1)));
				unsigned int cnt_2 = countbits_32(v_mask_left_2);
				v_mask_left_2 = ((v_mask_left_2 >> 1)) | v_mask_left_2;
				_MM_STOREU_MASK_INT((int *)(write_1), v_mask_left_2, v_i_16);
				write_1 += cnt_2; cnt_write_1 += cnt_2;

				SIMDMASKTYPE v_lt_weight_2 = _MM_CMP_LT_MASK(0xAAAA, _mm512_castsi512_ps(v_i_16), v_weight);
				SIMDMASKTYPE v_lt_id_2     = _MM_CMP_LT_MASK_INT(0x5555, v_i_16, v_id);
				SIMDMASKTYPE v_mask_right_2 = _mm512_kor(v_lt_weight_2, _mm512_kand(v_eq_weight_2, (v_lt_id_2<<1)));
				unsigned int cnt_right_2 = countbits_32(v_mask_right_2);
				v_mask_right_2 = ((v_mask_right_2 >> 1)) | v_mask_right_2;
				_MM_STOREU_MASK_INT((int *)(write_2), v_mask_right_2, v_i_16);
				write_2 += cnt_right_2; cnt_write_2 += cnt_right_2;
			}

			for(int i = tstart + num_aligned; i <= tend; i++)
			{
				if( verInd[i].weight > weight || (verInd[i].weight==weight && verInd[i].id>id))
				{
					*write_1 = verInd[i];
					write_1++; cnt_write_1++;
				}
				else if(verInd[i].weight < weight || (verInd[i].weight==weight && verInd[i].id<id))
				{
					*write_2 = verInd[i];
					write_2++; cnt_write_2++;
				}
			}

			// add the median at the end of write_1
			*write_1 = median;
                        write_1++; cnt_write_1++;

 			// copy back
			int cnt_write_2_aligned = cnt_write_2/8 * 8;

			for(int i = 0; i < cnt_write_2_aligned; i+=8)
			{
				_MM_STOREU_INT((int *)(write_1), _MM_LOAD_INT((int *)(neighbors + i))); write_1 += 8;
			}
			for(int i = cnt_write_2_aligned; i < cnt_write_2; i++)
			{
				*write_1 = neighbors[i];
				write_1++;
			}

			r = tstart + cnt_write_1 - 1;
#else
                        while ( p < r )
                        {
                            while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
                                p++;
                            while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
                                r--;

                            if(verInd[p].id==verInd[r].id)
                                p++;
                            else
                                if( p < r )
                                {
/*				    *(long long int *)(&temp) = *(long long int *)(verInd + p);*/
				    memcpy(&temp, verInd + r, sizeof(long long int));
				    *(long long int *)(verInd + p) = *(long long int *)(verInd + r);
/*				    *(long long int *)(verInd + r) = *(long long int *)(&temp);*/
				    memcpy(verInd + r, &temp, sizeof(long long int));
                                }
                        }
#endif
                        length=r-start+1;
                        if(length==k)
                            break;
                        else
                        {
                            if(length > k)
                                tend=r-1;
                            else
                                tstart=r+1;
                        }
                        if(tstart==tend)
                            break;

                    }

                    //nth_element(verInd+start,verInd+part,verInd+end,comparator);
                    sort(verInd+start,verInd+part,comparator);
                    //my_small_sort(verInd+start,(part - start));
                }
                break;

        default: sort(verInd+start,verInd+end,comparator);
    }
    return part;
}

int pivoting(Edge* verInd, int start, int end, int id, float weight, int* pindx,int step)
{
    int tstart,tend,r,p;

    tend=end-1;
    tstart=start;
    p=tstart;
    r=tend;
    while ( p < r )
    {
        while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
        {
            //if(p>tend)
                //break;
            p++;
        }
        while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
        {
            //if(r<tstart)
                //break;
            r--;
        }
        if(verInd[p].id==verInd[r].id)
            p++;
        else
            if( p < r )
            {
                /*temp.id= verInd[p].id;
                temp.weight=verInd[p].weight;

                verInd[p].id  = verInd[r].id;
                verInd[p].weight = verInd[r].weight;

                verInd[r].id = temp.id;
                verInd[r].weight = temp.weight;*/
                p++;
                r--;
            }
    }


    /*if(r<start)
        return start;
    else
    {
        if(r==tend)
        {
            sort(verInd+start,verInd+end,comparator);
            return end;
        }
        else
        {
            sort(verInd+start,verInd+r+1,comparator);
            return r+1;
        }
    }*/

    if(r<start)
        (*pindx)=start;
    else
    {
        if(r==tend)
            (*pindx)=end;
        else
            (*pindx)=r+1;
    }
    if(((*pindx)-start)>step)
        step=(*pindx)-start;
    return custom_sort(verInd,start,end,step,NULL,4,NULL);
}

int custom_sort(Edge* verInd, int start, int end, int step, int* order, int type, stack_pcl<Param>* st)
{
    int part=start+step;
    int tstart, tend,k,length;
    int id,p,r;
    float weight;
    Edge temp;

    switch(type)
    {
        case 1: break;
        case 2: part=end;
                sort(verInd+start,verInd+end,comparator);
                //csort(verInd,start,end);
                break;

        case 3: if(part>=end)
                {
                    part=end;
                    sort(verInd+start,verInd+end,comparator);
                }
                else
                {
                    int last=(end-start)/step;
                    //if(last>50)
                        //last=50;
                    //multiQselectIter(verInd,start,start,end-1,order,0,last-1,st);
                    multiQselect(verInd,start,start,end-1,order,0,last-1);
                    sort(verInd+start,verInd+part,comparator);
                }
                break;
        case 4: if(part>=end)
                {
                    part=end;
                    sort(verInd+start,verInd+end,comparator);
                }
                else
                {
                    tend=end-1;
                    tstart=start;
                    k=step+1;
                    while(true)
                    {
                        p=tstart;
                        r=tend;
                        weight = verInd[r].weight;
                        id=verInd[r].id;
                        while ( p < r )
                        {
                            while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
                                p++;
                            while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
                                r--;

                            if(verInd[p].id==verInd[r].id)
                                p++;
                            else
                                if( p < r )
                                {
                                    temp.id= verInd[p].id;
                                    temp.weight=verInd[p].weight;

                                    verInd[p].id  = verInd[r].id;
                                    verInd[p].weight = verInd[r].weight;

                                    verInd[r].id = temp.id;
                                    verInd[r].weight = temp.weight;
                                }
                        }

                        length=r-start+1;
                        if(length==k)
                            break;
                        else
                        {
                            if(length > k)
                                tend=r-1;
                            else
                                tstart=r+1;
                        }
                        if(tstart==tend)
                            break;
                    }

                    //nth_element(verInd+start,verInd+part,verInd+end,comparator);
                    sort(verInd+start,verInd+part,comparator);
                }
                break;

        default: sort(verInd+start,verInd+end,comparator);
    }
    return part;
}
