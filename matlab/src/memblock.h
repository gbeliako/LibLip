/**************************************************************************

    begin                : April 19 2004
	version				 : 1.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au

 * memblock.h -- service routines for memory allocation                    *
 * Used to take over OS heap, if many objects of a fixed size              *
 * need to be created, deleted and quickly accessed. This it to avoid      * 
 * OS keeping track of individual objects and associated overheads         *
 * Memblock implements a dynamic array to store all these objects, it's    *
 * own tracking system, and by having objects of the same size, reduces    *
 * overheads (by 10-100 times)                                             *
 *                                                                         *
 * An example of usage is to store a huge tree,  with nodes of equal size  *
 * Create a global variable                                                *
 *                                                                         *
 *  MemoryBlock<MyClass>	MB;                                          *
 *  UINT newcl, ref;  MyClass* cl;                                         *
 *  ref=MB.GetNextFree(num); equivalent of malloc(num*sizeof(Myclass))     *
 *	for(i=0;i<num;i++) {                                                 *  
 *		newcl=ref + i;         reference   to the i-th element         *
 *     	cl = MB.GetAt(newcl);  pointer to this class                   *
 *      // access members of cl                                            *
 *  }                                                                      *
 *  can free memory by                                                     *
 *  MB.FreeBlock(ref) ;                                                    *
 *  can release memory by calling ClearAll		                     *
 *                                                                         *
 *  © Gleb Beliakov, 2004								   *
 *                                                                         *
 * This program is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful, but     *
 * WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * General Public License for more details.                                *
 *                                                                         *
 * You should have received a copy of the GNU General Public License       *
 * along with this program; if not, write to the Free Software Foundation, *
 * Inc., 59 Temple Place Suite 330, Boston, MA 02111-1307 USA.             *
 ***************************************************************************/

#ifndef MEMORYBLOCK
#define MEMORYBLOCK

#include <cstdlib>
#include <malloc.h>
#include <memory.h>

#define UINT unsigned int

#define MB_IDX_SHF 20  // that's how we split the index into 2 parts: block and index within the block
#define MB_IDX_MASK ((1 << (MB_IDX_SHF)) - 1)
#define MB_BLK_SHF (32 - (MB_IDX_SHF))

/***  Macros for calculating the correct location of the node  ***/
#define MB_BLOCK(A) ((A) >> (MB_IDX_SHF))
#define MB_INDEX(A) ((A) & (MB_IDX_MASK))
#define MB_INDEXB(A,B) (((A) << (MB_IDX_SHF)) + B)

/***  The upper limits - 4GB  ***/
#define MB_MAX_NODES 0xFFFFFFFF

/***  Define the ranges  ***/
#define MB_MAX_INDEXES (1 << (MB_IDX_SHF))
#define MB_MAX_BLOCKS (((MB_MAX_NODES) / (MB_MAX_INDEXES) + 1))

#define MB_BADINDEX		0xFFFFFFFF
#define MB_SPECIALINDEX 0xEFFFFFFF
#define MB_BLOCKSIZE	0x7FFFF

#define	HasFree(B) ( !((B) & 0x1) )
#define	HasAnyFree(B) ( (B!=0xFFFFFFFE) )
#define SetFree(B) ( ((B) &= 0xFFFFFFFE) )

#define SetFreeI(B,r)   { (B) &= (~(0x1 << r)) ;(B) &= 0xFFFFFFFE; } 


inline int  SetOccupied(UINT &B, short i) { 
	B |= (0x1 << i); 
	if(HasAnyFree(B)) {SetFree(B); return 0; } 
	else  {B |= 0x00000001; return 1;}
};

inline short WhichFree(UINT B) {
	for(short   i=1;i<32;i++)
		if(!((B>>i) & 0x1) ) return i;
	return 0;
};

template <class T>
class MemBlock {
public:

    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;


	value_type * m_data;
//	UINT m_index[1024];
	UINT	m_NextAvail,m_temp;

	short i,j,k;
	//0 means free block, 1 means occupied

	MemBlock() {       // 7FFF= 31*32*32 +31*32 +31
		m_data=(value_type*) calloc(MB_BLOCKSIZE,sizeof(T) ); // 31000 blocks of size T, >64kb
//		for(short  i=0;i<1024;i++) m_index[i]=0;
		memset(m_data,0xFF,MB_BLOCKSIZE*sizeof(T));
		m_NextAvail=0;
	};

	~MemBlock() {
		free(m_data);
	};

	UINT	GetNextFree() {
		if(m_NextAvail>=MB_BLOCKSIZE) return MB_BADINDEX;
		m_temp=m_NextAvail;
		m_NextAvail++;
		return m_temp;
	};

	UINT	GetNextFree(int M) {
		if(m_NextAvail+M>=MB_BLOCKSIZE) return MB_BADINDEX;
		m_temp=m_NextAvail;
		m_NextAvail+=M;
		return m_temp;
	};


	inline int	IsFree() { return (m_NextAvail<MB_BLOCKSIZE); };
	inline int	IsFreeM(int M) { return (m_NextAvail+M < MB_BLOCKSIZE); };
	inline void	FreeBlock(UINT B) { memset(m_data+B,0xFF,sizeof(T)); };

	T* GetAt(UINT B) { return (T*) (m_data+B); };

	void SetAt(UINT B, T* Value) { memcpy(m_data + B, Value, sizeof(T)); };

	inline reference operator()(UINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (UINT B) const { *((T*)(GetAt(B))); };

	void ClearAll() {	m_NextAvail=0;};
};




template <class T>
class MemoryBlock 
{
public:
    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;

  MemBlock<T>** block;

  UINT nodeCount, emptyBlocks, currentBlock;
  int valid;

	MemoryBlock(void)
	{
		block =  (MemBlock<T>** ) calloc(MB_MAX_BLOCKS, sizeof(MemBlock<T>*));
//		assert(block != NULL);
		nodeCount = emptyBlocks = 0;
/***  The starting point is 0 but set to -1 because the
      _createNextBlock will increament the value
      before using it
***/
		currentBlock = (UINT) -1;
		_createNextBlock();
		valid=1;
	};

	~MemoryBlock(void)
	{
		for(UINT  loop = currentBlock + emptyBlocks; loop > 0; loop--) delete (block[loop]);
/***  To free the first block!  ***/
	  delete(block[0]);
	  free(block);
	  valid=0;
	};

	UINT	GetNextFree() {
		nodeCount++;
		UINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree());
				return loop;
			}
		// no space left
		if(currentBlock < MB_MAX_BLOCKS-2) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree());
			return loop;
		}
		nodeCount--;
		//exit(20);
		return MB_BADINDEX;
	};

	UINT	GetNextFree(int M) {
		nodeCount+=M;
		UINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFreeM(M)) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree(M));
				return loop;
			}
		// no space left
		if(currentBlock < MB_MAX_BLOCKS-2) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree(M));
			return loop;
		}
		nodeCount--;
		//exit(20);
		return MB_BADINDEX;
	};

	inline void	FreeBlock(UINT B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
	//	B=MB_BADINDEX;
	};

	inline int	IsFree() { 
		if(currentBlock < MB_MAX_BLOCKS-1) return 1;
		for(UINT  loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) return 1;
		return 0; 
	};

	inline T* GetAt(UINT B) { return block[MB_BLOCK(B)]->GetAt(MB_INDEX(B)); };

	inline void SetAt(UINT B, T* Value) { block[MB_BLOCK(B)]->SetAt(MB_INDEX(B),Value); };

	inline reference operator()(UINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (UINT B) const { *((T*)(GetAt(B))); };

	inline void _createNextBlock()
	{
	  currentBlock++;

	  if(emptyBlocks == 0)
	  {
		block[currentBlock] = new MemBlock<T>;
//		assert(block[currentBlock] != NULL);
	  }
	  else
		emptyBlocks--;
	}

	void ClearAll() 
	{
//		cout << "commited blocks " <<currentBlock + emptyBlocks<<" of size "<<block[0]->m_NextAvail <<endl;

		for(UINT  loop = currentBlock + emptyBlocks; loop > 0; loop--) block[loop]->ClearAll();
/***  To free the first block!  ***/
	  block[0]->ClearAll();
	};

	inline int IsValid() {return valid;}
};




/*-------------old version for lists------------------------*/
template <class T>
class MemBlockE {
public:

    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;


	value_type * m_data;
	UINT m_index[1024];

	short i,j,k;
	//0 means free block, 1 means occupied

	MemBlockE() {       // 7FFF= 31*32*32 +31*32 +31
		m_data=(value_type*) calloc(0x7FFF,sizeof(T) ); // 31000 blocks of size T, >64kb
		for(short   i=0;i<1024;i++) m_index[i]=0;
		memset(m_data,0xFF,0x7FFF*sizeof(T));
	};

	~MemBlockE() {
		free(m_data);
	};

	UINT	GetNextFree() {
		//short   i,j,k;
		i=WhichFree(m_index[0]);
		j=WhichFree(m_index[i]);
		k=WhichFree(m_index[i*32+j]);
		if(SetOccupied(m_index[i*32+j],k))
			if(SetOccupied(m_index[i],j))
				SetOccupied(m_index[0],i);
		return GetAddress();// i,j,k);
	};

	void	FreeBlock(UINT B)	{
		//short i,j,k;
		GetIJK(B);//,i,j,k);
		SetFreeI(m_index[i*32+j],k); 
		SetFreeI(m_index[i],j); 
		SetFreeI(m_index[0],i);
	};

	UINT	GetAddress(){ //short i, short j, short k)	{
		UINT r=(i-1);
		r *= 1024;
		r = r+ (j-1)*32 + k-1;
		//return ((i-1)*32*32+(j-1)*32 +k-1 ); //sizeof(MyStruct_t)*
		return r;
	};

	void	GetIJK(UINT B)//, short& i, short& j, short &k)
	{
		div_t t=div(B,32); ///sizeof(MyStruct_t)
		j=t.quot; k=t.rem+1;
		t=div(j,32);
		j=t.rem+1; i=t.quot+1;
	};

	inline int	IsFree() { return HasFree(m_index[0]); };

	T* GetAt(UINT B) { return (T*) (m_data+B); };

	void SetAt(UINT B, T* Value) { memcpy(m_data + B, Value, sizeof(T)); };

	inline reference operator()(UINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (UINT B) const { *((T*)(GetAt(B))); };

	void ClearAll() {
		for(short   i=0;i<1024;i++) m_index[i]=0;
	};
};




template <class T>
class MemoryBlockE 
{
public:
    typedef         T   value_type;
    typedef         T*  pointer;
    typedef         T&  reference;
    typedef const   T&  const_reference;

  MemBlockE<T>** block;

  UINT nodeCount, emptyBlocks, currentBlock;
  int valid;

	MemoryBlockE(void)
	{
		block =  (MemBlockE<T>** ) calloc(MB_MAX_BLOCKS, sizeof(MemBlockE<T>*));
//		assert(block != NULL);
		nodeCount = emptyBlocks = 0;
/***  The starting point is 0 but set to -1 because the
      _createNextBlock will increament the value
      before using it
***/
		currentBlock = (UINT) -1;
		_createNextBlock();
		valid=1;
	};

	~MemoryBlockE(void)
	{
		for(UINT   loop = currentBlock + emptyBlocks; loop > 0; loop--) delete (block[loop]);
/***  To free the first block!  ***/
	  delete(block[0]);
	  free(block);
	  valid=0;
	};

	UINT	GetNextFree() {
		nodeCount++;
		UINT loop;
		for(loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) {
				loop = MB_INDEXB(loop, block[loop]->GetNextFree());
				return loop;
			}
		// no space left
		if(currentBlock < MB_MAX_BLOCKS-2) {
			_createNextBlock();
			loop = MB_INDEXB(loop, block[currentBlock]->GetNextFree());
			return loop;
		}
		nodeCount--;
		exit(20);
		return MB_BADINDEX;
	};

	inline void	FreeBlock(UINT& B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
		B=MB_BADINDEX;
		//if(BLOCK(B) == currentBlock  && 
	};
	inline void	FreeBlockC(UINT B)	{
		nodeCount--;
		block[MB_BLOCK(B)]->FreeBlock(MB_INDEX(B));
	};


	inline int	IsFree() { 
		if(currentBlock < MB_MAX_BLOCKS-1) return 1;
		for(UINT   loop=0; loop <= currentBlock; loop++)
			if(block[loop]->IsFree()) return 1;
		return 0; 
	};

	inline T* GetAt(UINT B) { return block[MB_BLOCK(B)]->GetAt(MB_INDEX(B)); };

	inline void SetAt(UINT B, T* Value) { block[MB_BLOCK(B)]->SetAt(MB_INDEX(B),Value); };

	inline reference operator()(UINT B) { return *((T*)(GetAt(B)));	};
    inline const_reference operator() (UINT B) const { *((T*)(GetAt(B))); };

	inline void _createNextBlock()
	{
	  currentBlock++;

	  if(emptyBlocks == 0)
	  {
		block[currentBlock] = new MemBlockE<T>;
//		assert(block[currentBlock] != NULL);
	  }
	  else
		emptyBlocks--;
	}

	void ClearAll() 
	{
	  for(UINT   loop = currentBlock + emptyBlocks; loop > 0; loop--) block[loop]->ClearAll();
/***  To free the first block!  ***/
	  block[0]->ClearAll();
	};

	inline int IsValid() {return valid;}
};



#endif
