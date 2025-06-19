/*
 *  MatlabFileHandle.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 7/26/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#ifndef _MATLABFILEHANDLE_H_
#define _MATLABFILEHANDLE_H_

#include <Main.H>
#include <Vector.H>
#include <Matrix.H>
#include <FileHandle.H>

/* Matlab header */
/* FORMAT = 1000 on PPC and FORMAT = 0000 on Intel */
#define FORMAT 0000

/* A struct that takes the Matlab header */

class MatlabHeader
	{
	public:
		
		MatlabHeader(long type, long rows, long cols, long imagf, long namelen) 
		{
			m_spec[0] = type;
			m_spec[1] = rows;
			m_spec[2] = cols;
			m_spec[3] = imagf;
			m_spec[4] = namelen;
		}
		
		long m_spec[5]; /* type, rows, cols, imagf, namelen */
	};


/* A class of MatlabFileHandle */

class MatlabFileHandle: public FileHandle 
	{
	public:
		
		/* Forbid default constructor */
		MatlabFileHandle();
		
		/* Constructor with file name and opening mode */
		MatlabFileHandle(const std::string &a_fileName, openingModes a_openingMode) :
		FileHandle::FileHandle(a_fileName) 
		{
			if (a_openingMode==OPEN_RD)
				m_fstream.open(m_fileName.data(), std::ios::binary | std::ios::in);
			else if (a_openingMode==OPEN_WR)
				m_fstream.open(m_fileName.data(), std::ios::binary | std::ios::out);
		}
		
		/* Write a scalar */
		void write(const double a_value, const std::string &a_name) 
		{
			MatlabHeader header(FORMAT,1,1,0,a_name.size() + 1);
			
			m_fstream.write((char*)header.m_spec, sizeof(MatlabHeader));
			m_fstream.write((char*)a_name.data(), header.m_spec[4]);
			m_fstream.write((char*)&a_value, sizeof(double));
			m_fstream.flush();
		}

		/* Write a vector */
		void write(const Vector<double> &a_vector, const std::string &a_name) 
		{	
			int nlen = a_vector.length();
			MatlabHeader header(FORMAT,nlen,1,0,a_name.size() + 1);
			
			m_fstream.write((char*)header.m_spec, sizeof(MatlabHeader));
			m_fstream.write((char*)a_name.data(), header.m_spec[4]);
			m_fstream.write((char*)a_vector.getPointer(), sizeof(double)*nlen);
			m_fstream.flush();	
		}

		/* Write a Matrix */
		void write(const Matrix<double> &a_matrix, const std::string &a_name) 
		{	
			int nrows = a_matrix.rows();
			int ncols = a_matrix.cols();
			MatlabHeader header(FORMAT,nrows,ncols,0,a_name.size() + 1);
			
			m_fstream.write((char*)header.m_spec, sizeof(MatlabHeader));
			m_fstream.write((char*)a_name.data(), header.m_spec[4]);
			m_fstream.write((char*)a_matrix.getPointer(), sizeof(double)*nrows*ncols);
			m_fstream.flush();
		}
		
	};




#endif
