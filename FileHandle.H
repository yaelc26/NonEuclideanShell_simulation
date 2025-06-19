/*
 *  FileHandle.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 7/26/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

/* 
	A FileHandle is an abstract class for a file handle.
	It has a filename associated with an fstream. It really is a base class.
	We derive from it MatlabFileHandle and TextFileHandle.
*/

#ifndef _FILEHANDLE_H_
#define _FILEHANDLE_H_

#include "Main.H"
#include <string>
#include "Errors.H"

class FileHandle
		{
		public:
			/* Opening modes */
			enum openingModes {OPEN_RD, OPEN_WR};
			
			/* Default constructor: not allowed */
			FileHandle();
			
			/* Constructor with file name and opening mode */
			FileHandle(const std::string &a_fileName) :
				m_fileName(a_fileName)
			{}
			
			/* Copy constructor: not allowed */
			FileHandle(const FileHandle &);
			
			/* Close the stream */
			void close() {m_fstream.close();}
	
		protected:

			std::string    m_fileName;
			std::fstream   m_fstream; 

		};



#endif
