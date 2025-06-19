/*
 *  TextFileHandle.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/13/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#ifndef _TEXTFILEHANDLE_H_
#define _TEXTFILEHANDLE_H_

#include <Main.H>
#include <FileHandle.H>



class TextFileHandle: public FileHandle 
	{
	public:
		
		/* Forbid default constructor */
		TextFileHandle();
		
		/* Constructor with file name and opening mode */
		TextFileHandle(const std::string &a_fileName, openingModes a_openingMode) :
		FileHandle::FileHandle(a_fileName) 
		{
			if (a_openingMode==OPEN_RD)
				m_fstream.open(m_fileName.data(), std::ios::in);
			else if (a_openingMode==OPEN_WR)
				m_fstream.open(m_fileName.data(), std::ios::out);
		}
		
		
		/* Write a double floats and ints */
		void write(char*  a_ptr)    {m_fstream << a_ptr;}
		void write(int    a_value)  {m_fstream << a_value;}
		void write(double a_value)  {m_fstream << a_value;}

		void newline()             {m_fstream << std::endl;}
		void tab()                 {m_fstream << "\t";}

		/* Read doubles, floats and ints */
		void read(int&    a_ref) {m_fstream >> a_ref;}
		void read(double& a_ref) {m_fstream >> a_ref;}
		void read(float&  a_ref) {m_fstream >> a_ref;}
	};

#endif

