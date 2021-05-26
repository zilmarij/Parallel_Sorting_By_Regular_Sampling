
#include<iostream>
#include<math.h>
#include<cstring>
#include "mpi.h"
#include "sort.h"



//using namespace std;
int count;
int size;  int leng;
int rank; int rem;
pSort::dataType* f_buf;
int offset;



void merge(pSort::dataType* result, pSort::dataType* left, pSort::dataType* right, int leftLen, int rightLen);
pSort::dataType* slice(pSort::dataType* rbuf, int start, int end);
void MERGe(pSort::dataType* rbuf, int len);

void place(pSort::dataType* rbuf, int place_val, int c);
void RADIx(pSort::dataType* rbuf, int c);

int part(pSort::dataType* rbuf, int low, int high);
void QUICk(pSort::dataType* rbuf, int l, int h);


void QUICK2(pSort::dataType* rbuf);			//local pivot extraction, sharing with P0; obtaining pivots from P0 
void QUICK3(int* recvd_indices, pSort::dataType*);   //partitioniing using rxed pivots from P0

void QUICKint(int* buff, int l, int h);
int partint(int* buf, int low, int high);


void pSort::init()
{
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

void pSort::sort(pSort::dataType* rbuf, int l, pSort::SortType type)
{
	//std::cout << " leng is " << l;
	leng = l * size;
	//count = (int)(leng / size);			//count of array elements to be sent to each process
	//rem = (int)(leng % size);			//number of processes to send count+1 elements
	rem = 0; count = l;


	
	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{

			if (type == 1 || type == 0)
			{
				QUICk(rbuf, 0, count-1);
			}
			else if (type == 2)
			{
				MERGe(rbuf, count);
			}
			else if (type == 3)
			{
				RADIx(rbuf, count);
			}

		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			QUICK2(rbuf);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 100)
	{
		std::cout << "Precs1 \n";
		for (int i = 915; i < 930; i++)
		{
			std::cout << rbuf[i].key << "\n";
			
		}
	}
}



void pSort::close()
{
	MPI_Finalize();
	return;
}


void INSERTION(int* rbuf, int c)
{
	for (int i = 1; i < c; i++)
	{
		int tmp = rbuf[i], j;

		for (j = i; j >= 1 && tmp < rbuf[j - 1]; j--)
			rbuf[j] = rbuf[j - 1];

		rbuf[j] = tmp;
	}
}


void place(pSort::dataType* rbuf, int place_val, int c)
{
	for (int z = 0; z < size; z++)
	{
		if (rank == z)
		{
			pSort::dataType* temp = new pSort::dataType[c + 1];
			int* tmp = new int[10];  temp[0].key = 0; int neg = 0;
			int ng = 0;
			for (int i = 0; i < 10; i++)  //arr of place values
			{
				tmp[i] = 0;
			}
			for (int i = 0; i < c; i++)
			{
				temp[i].key = i;
			}

			for (int i = 0; i < c; i++)
			{
				if (rbuf[i].key < 0)
				{
					neg++;
				}
				else
				{
					int v = (rbuf[i].key / place_val) % 10;
					tmp[(rbuf[i].key / place_val) % 10] += 1;  //record number of occurences of digit i
				}
				
			}
			tmp[0] += neg; ng = neg;
			for (int i = 1; i < 10; i++)
			{
				tmp[i] += tmp[i - 1]; //displacement of digit i
			}

			for (int i = c - 1; i >= 0; i--)
			{
				if (rbuf[i].key < 0)
				{
					int x = neg - 1; neg--;
					std::memcpy(temp + x, rbuf + i, sizeof(rbuf));
				}
				else
				{
					int x = tmp[(rbuf[i].key / place_val) % 10] - 1;
					std::memcpy(temp + x, rbuf + i, sizeof(rbuf));

					
					tmp[(rbuf[i].key / place_val) % 10] -= 1;
				}
			}

			for (int i = 0; i < c; i++)
			{
				
				std::memcpy(rbuf + i, temp + i, sizeof(rbuf));
			}
			if (place_val == 1000000000)
			{
				QUICk(rbuf, 0, ng);
			}

		}
	}
}


void RADIx(pSort::dataType* rbuf, int c)
{
	for (int i = 0; i < 10; i++)		//10 because the range of 32-bit int allows these many digits
	{
		place(rbuf, pow(10, i), c);
	}
}


pSort::dataType* slice(pSort::dataType* rbuf, int start, int end)
{
	pSort::dataType* result = new pSort::dataType[(end - start)];
	int i;
	for (i = start; i < end; i++)
	{
		result[i - start] = rbuf[i];
	}

	return result;
}

void merge(pSort::dataType* result, pSort::dataType* left, pSort::dataType* right, int leftLen, int rightLen)
{
	int i = 0, j = 0;
	while (i < leftLen && j < rightLen)
	{
		if (left[i].key < right[j].key)
		{
			result[i + j] = left[i];
			i++;
		}
		else
		{
			result[i + j] = right[j];
			j++;
		}
	}

	for (; i < leftLen; i++)
	{
		result[i + j] = left[i];
	}
	for (; j < rightLen; j++)
	{
		result[i + j] = right[j];
	}


}

void MERGe(pSort::dataType* rbuf, int len)
{
	if (len <= 1)
	{
		return;
	}
	pSort::dataType* left = slice(rbuf, 0, len / 2 + 1);
	pSort::dataType* right = slice(rbuf, len / 2, len);

	MERGe(left, len / 2);
	MERGe(right, len - (len / 2));

	merge(rbuf, left, right, len / 2, len - (len / 2));
}



int part(pSort::dataType* rbuf, int low, int high)
{
	int pivot = rbuf[high].key;    // pivot 
	int i = (low - 1);  // Index of smaller element 

	for (int j = low; j <= (high - 1); j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (rbuf[j].key <= pivot)
		{
			i++;
			pSort::dataType t = rbuf[i];
			rbuf[i] = rbuf[j];
			rbuf[j] = t;
		}
	}
	pSort::dataType t = rbuf[i + 1];
	rbuf[i + 1] = rbuf[high];
	rbuf[high] = t;
	return (i + 1);
}

void QUICk(pSort::dataType* rbuf, int l, int h)
{

	if (l < h)
	{
		int p = part(rbuf, l, h);

		QUICk(rbuf, l, p - 1);
		QUICk(rbuf, p + 1, h);

	}
}

void QUICKint(int* buff, int l, int h)
{
	if (l < h)
	{
		int p = partint(buff, l, h);

		QUICKint(buff, l, p - 1);
		QUICKint(buff, p + 1, h);

	}
}
int partint(int* buf, int low, int high)
{
	int pivot = buf[high];    // pivot 
	int i = (low - 1);  // Index of smaller element 

	for (int j = low; j <= (high - 1); j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (buf[j] <= pivot)
		{
			i++;
			int t = buf[i];
			buf[i] = buf[j];
			buf[j] = t;
		}
	}
	int t = buf[i + 1];
	buf[i + 1] = buf[high];
	buf[high] = t;
	return (i + 1);
}

void QUICK2(pSort::dataType* rbuf)
{
	int* local_indices = new int[size];

	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			//std::cout << " quick2 for " << rank;
			for (int z = 0; z < size; z++)
			{
				int p = size * size;
				int y = floor(z * leng / p);
				local_indices[z] = rbuf[y].key;
			}
		}
	}

	int* pivots_recvd = new int[(size * size)]; //for P0 to receive from others

	int* recvd_indices = new int[(size - 1)];    //for recving from P0

	if (rank == 0)
	{
		for (int i = 0; i < size; i++)
		{
			pivots_recvd[i] = local_indices[i];
			
		}
	}


	for (int i = 1; i < size; i++)
	{
		if (rank == i)		//processes send their pivots to P0
		{

			MPI_Send(local_indices, size, MPI_INT, 0, 2, MPI_COMM_WORLD);

		}
	}

	if (rank == 0)
	{
		for (int i = 1; i < size; i++)		//P0 recvs pivots frm other processes
		{
			MPI_Recv(pivots_recvd + (i * size), size, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}

		int  l = 0; int h = (size * size) - 1;

		if (l < h)		//P0 sorts recvd pivots
		{
			int p = partint(pivots_recvd, l, h);

			QUICKint(pivots_recvd, l, p - 1);
			QUICKint(pivots_recvd, p + 1, h);
		}


		int* pivots_sent = new int[(size - 1)];
		int t = (floor(size / 2)) - 1;

		for (int i = 1; i < size; i++)
		{
			int y = i * size + t;
			pivots_sent[i - 1] = pivots_recvd[y];
		}
		
		for (int i = 0; i < size - 1; i++)
		{
			std::memcpy(recvd_indices + i, pivots_sent + i, sizeof(recvd_indices));
		}

		for (int i = 1; i < size; i++)
		{
			MPI_Send(pivots_sent, size - 1, MPI_INT, i, 3, MPI_COMM_WORLD);
		}

		
	}

	for (int i = 1; i < size; i++)
	{
		if (rank == i)
		{
			
			MPI_Recv(recvd_indices, size - 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	

	QUICK3(recvd_indices, rbuf);


}

void QUICK3(int* recvd_indices, pSort::dataType* rbuf)
{


	MPI_Datatype Type;
	
	f_buf = (pSort::dataType*)calloc(10, sizeof(pSort::dataType));   //final buffer, to hold partitions
	int blocklengs[2]; MPI_Aint indices[2]; MPI_Datatype old_types[2];
	blocklengs[0] = 1; blocklengs[1] = LOADSIZE;
	old_types[0] = MPI_INT; old_types[1] = MPI_CHAR;
	MPI_Get_address(&((*f_buf).key), &indices[0]);
	MPI_Get_address(&((*f_buf).payload), &indices[1]);
	indices[1] = indices[1] - indices[0];
	indices[0] = 0;
	MPI_Type_create_struct(2, blocklengs, indices, old_types, &Type);
	MPI_Type_commit(&Type);

	offset = 0; int add = 0;
	MPI_Status status;
	int send_count = 0;


	
	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			int c = 0;
			int begin = 0;
			for (int j = 0; j < size; j++)	 //partition and send
			{
				send_count = 0;
				pSort::dataType* inter; 
				begin = c;

				//if ((i < rem && c > count) || (i >= rem && c > (count - 1)))	//last part not found, end of items
				if(c==count)
				{
					send_count = 0;
				}

				else if (j == size - 1)		//exhausted the pivots which are one less than the processes, made it to the last partitioin
				{
					/*if (i < rem)
					{
						send_count = count + 1 - c;
						//std::cout << " \n extended < rem" << send_count;

					}
					else
					{*/
					send_count = count - c; //c = count;
						//std::cout << " \n extended " << send_count;
					//}
				}
				else
				{
					while (true)
					{
						//if ((i < rem && c == count) || (i >= rem && c == (count - 1)))  //end of items while checking
						if(c==count)
						{
							break;
						}
						if (rbuf[c].key <= recvd_indices[j])	//c will determine the partition count
						{
							//std::cout << "\n found for it " <<j<< " " << rbuf[c].key;
							c++;
							send_count++;

						}
						else
						{
							//std::cout << "\n send count for j " <<j<< " " << send_count << std::endl;
							break;
						}
					}
				}
				inter = new pSort::dataType[send_count+1];
				inter[0].key = 0;
				if (send_count == 0) { inter = nullptr; }

				for (int g = 0; g < send_count; g++)		//copy the partition to send fwd
				{
					inter[g].key = rbuf[begin + g].key;
					for (int k = 0; k < LOADSIZE; k++)
					{
						inter[g].payload[k] = rbuf[begin + g].payload[k];
					}
					//memcpy(inter + g, rbuf + begin + g, sizeof(inter));
					//std::cout << " memcpy " << inter[g].key << " "<< " rank "<< j;
				}

				if (rank == j)
				{
					if (inter != nullptr)
					{
						pSort::dataType* np;
						if (np = (pSort::dataType*)realloc(f_buf, (offset + send_count + 1) * sizeof(pSort::dataType)))
						{
							f_buf = np;
							//std::cout << " realloc succsful for" << rank;
						}
						for (int r = 0; r < send_count; r++)   //to self
						{
							f_buf[r + offset].key = rbuf[begin + r].key;
							for (int k = 0; k < LOADSIZE; k++)
							{
								f_buf[r + offset].payload[k] = rbuf[begin + r].payload[k];
							}
							//std::memcpy(f_buf + r + offset, rbuf + begin + r, sizeof(f_buf));

						}
						offset += send_count;
					}
				}
				else  //senders' rank != receivers' rank
				{
					MPI_Send(inter, send_count, Type, j, i, MPI_COMM_WORLD);
				}
				//delete[] inter;
			}

		}

		else if (rank != i)
		{
			MPI_Probe(i, i, MPI_COMM_WORLD, &status);			
			MPI_Get_count(&status, Type, &add);
			//std::cout << rank << " recving " << add << " frm " << i;
			
			pSort::dataType* np;
		    if (np = (pSort::dataType*)realloc(f_buf, (offset+add+1) * sizeof(pSort::dataType)))
			{
				f_buf = np;
				//std::cout << " ReallocC succsful for " << rank;
			}
			//f_buf = (pSort::dataType*)realloc(f_buf, (offset) * sizeof(pSort::dataType));
			
			MPI_Recv(f_buf + offset, add, Type, i, i, MPI_COMM_WORLD, &status);
			//std::cout << "received " << add << " indices of data by " << rank <<" from " << i <<" offset + add " << offset+add<< "\n";
			offset += add;

		}

	}

	
	
	for (int i = 0; i < size; i++)
	{
		if (rank == i)
		{
			QUICk(f_buf, 0, offset - 1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	int ctr = count;  //index frm whr to send fwd, pertains to f_buf
	int diff = offset - count;    //how much to send fwd or receive
	int dif[1];
	dif[0] = diff;
	int recv_diff[1]; // = new int[1];
	int flag = 0; int flag2 = 0;
	int v;
	if (offset <= count)
	{
		v = offset;
	}
	else if (offset > count)
	{
		v = count;
	}

	if (rank == 0)  //if diff=0, no sending or receiving reqd
	{
		for (int i = 0; i < v; i++)
		{
			rbuf[i].key = f_buf[i].key;
			for (int ii = 0; ii < LOADSIZE; ii++)
			{
				rbuf[i].payload[ii] = f_buf[i].payload[ii];
			}
			
		}
		flag2 = v; //to send data from f_buf, v to get data into rbuf
		
	}
	for (int i = 0; i < size - 1; i++)
	{
		if (rank == i)
		{
			
			dif[0] = diff;
			MPI_Send(dif, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);   // tell the nxt p abt excess or deficit

			if (diff > 0)   //send fwd the excessive data
			{
				MPI_Send(f_buf + flag2, diff, Type, rank + 1, rank, MPI_COMM_WORLD);  //rbuf is set in this case
			}
			else if (diff < 0)  //get the defict data
			{
				int gt = count - v;
				MPI_Recv(rbuf + v, gt, Type, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
				

			}

		}

		else if (rank == i + 1)
		{

			MPI_Status status;
			MPI_Recv(recv_diff, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
			
			if (recv_diff[0] > 0) //prev process has excessive data, recv it
			{
				MPI_Status status;
				MPI_Recv(rbuf, count, Type, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, Type, &add);
				

				diff = count - add;  //how much more to add from f_buf
				ctr = 0;

				if (diff < offset) { v = diff; }  //if offset > remaining reqd data
				else if (diff >= offset) { v = offset; }
				for (int k = 0; k < v; k++)
				{
					rbuf[add + k].key = f_buf[k].key;
					
					for (int ii = 0; ii < LOADSIZE; ii++)
					{
						rbuf[add + k].payload[ii] = f_buf[k].payload[ii];
					}
					ctr++;
				}
				diff = offset - (count - add);  //how much to get or give
				flag2 = v;  //offset for f_buf to send across
				v = add + v;  //offset for rbuf to recv
			}
			else if (recv_diff[0] < 0)  //the prev process is deficit of data, send it
			{
				int rd = recv_diff[0];
				rd = rd * (-1); flag = rd;
				
				MPI_Send(f_buf, rd, Type, rank - 1, rank, MPI_COMM_WORLD);
				//v = offset;
				if (count ==(offset-rd))    //remaining elements will suffice its need
				{
					v = offset;
				}
				else if (count < (offset - rd))
				{
					v = count+rd;
				}
				else if (count > (offset - rd))
				{
					v = offset;
				}
				int abc = 0;
				for (int k = rd; k < v; k++)
				{
					rbuf[abc].key = f_buf[k].key;
					for (int ii = 0; ii < LOADSIZE; ii++)
					{
						rbuf[abc].payload[ii] = f_buf[k].payload[ii];
					}
					//memcpy(rbuf + k + add, f_buf + k, sizeof(rbuf));
					ctr++; abc++;
				}
				diff = offset - rd - count; 
				flag2 = v; //whr to fwd data from f_buf
				v = abc;    //whr to add data in rbuf

				
			}
			else if (recv_diff[0] == 0)
			{
				
				for (int k = 0; k < v; k++)
				{
					rbuf[k].key = f_buf[k].key;
					
					for (int ii = 0; ii < LOADSIZE; ii++)
					{
						rbuf[k].payload[ii] = f_buf[k].payload[ii];
					}
				}
				diff = offset - count; 
				flag2 = v; //to send data from f_buf, v to add data to rbuf

			}
		}

	}


}
