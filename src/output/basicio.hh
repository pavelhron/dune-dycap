#ifndef DUNE_DYCAP_BASICIO_HH
#define DUNE_DYCAP_BASICIO_HH
#include <fstream>
#include <numeric>
#include <inttypes.h>
#include <stdint.h>
#include <iostream>

/** \brief Writing of solution vector.
 */
template< typename GFS, typename U >
void writeGFSVector(std::ofstream & file, const GFS & gfs, const U & u)
{
  typedef uint64_t IOSizeType;
  const IOSizeType size = gfs.size();

  file.write((char*)&size,sizeof(IOSizeType));
  for (auto it = u.begin(); it!=u.end();it++)
    file.write((char*)&(*it),sizeof(typename U::field_type));
}

/** \brief Reading of stored solution vector.
 */
template< typename GFS, typename U >
void readGFSVector
(std::ifstream & file, const GFS & gfs, U & u)
{
  typedef uint64_t IOSizeType;
  const IOSizeType gfs_size = gfs.size();

  IOSizeType size;
  file.read((char*)&size,sizeof(IOSizeType));

  if(size != gfs_size){
    std::cerr << "Number of DOFs indicated in file does not match number "
      "of DOFs in grid function space!" << std::endl;

    u.base().resize(size);
  }

  for (auto it = u.begin(); it!=u.end();it++)
    file.read((char*)&(*it),sizeof(typename U::field_type));
}

/** \brief Writing of std solution vector.
 */
template< typename V >
void writeSTDVector(std::ofstream & file, const V & v)
{
  typedef uint64_t IOSizeType;
  const IOSizeType size = v.size();
  file.write((char*)&size,sizeof(IOSizeType));

  typedef typename V::value_type VT;
  for(IOSizeType i=0; i<size; ++i)
    file.write((char*)&(v[i]),sizeof(VT));
}

/** \brief Reading of stored std solution vector.
 */
template< typename V >
void readSTDVector(std::ifstream & file, V & v)
{
  typedef uint64_t IOSizeType;
  IOSizeType size;
  file.read((char*)&size,sizeof(IOSizeType));
  v.resize(size);

  typedef typename V::value_type VT;
  for(IOSizeType i=0; i<size; ++i)
    file.read((char*)&(v[i]),sizeof(VT));
}

#endif
