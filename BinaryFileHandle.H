/*
 *  BinaryFileHandle.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/11/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#ifndef _BINARYFILEHANDLE_H_
#define _BINARYFILEHANDLE_H_

#include <Main.H>
#include <FileHandle.H>



class BinaryFileHandle: public FileHandle 
	{
	public:
		
		/* Forbid default constructor */
		BinaryFileHandle();
		
		/* Constructor with file name and opening mode */
		BinaryFileHandle(const std::string &a_fileName, openingModes a_openingMode) :
		FileHandle::FileHandle(a_fileName) 
		{
			if (a_openingMode==OPEN_RD)
				m_fstream.open(m_fileName.data(), std::ios::binary | std::ios::in);
			else if (a_openingMode==OPEN_WR)
				m_fstream.open(m_fileName.data(), std::ios::binary | std::ios::out);
		}
		
		
		/* Write a double floats and ints */
		void write(double* a_ptr, int a_length=1) {m_fstream.write((char*)a_ptr, sizeof(double)*a_length);}
		void write(float*  a_ptr, int a_length=1) {m_fstream.write((char*)a_ptr, sizeof(float)*a_length);}
		void write(int*    a_ptr, int a_length=1) {m_fstream.write((char*)a_ptr, sizeof(int)*a_length);}
		
		/* read */
		void read(double* a_ptr, int a_length=1)  {m_fstream.read((char*)a_ptr, sizeof(double)*a_length);}
		void read(float*  a_ptr, int a_length=1)  {m_fstream.read((char*)a_ptr, sizeof(float)*a_length);}
		void read(int*    a_ptr, int a_length=1)  {m_fstream.read((char*)a_ptr, sizeof(int)*a_length);}
		
	};

#endif

