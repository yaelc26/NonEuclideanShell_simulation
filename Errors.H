/*
 *  Errors.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 8/12/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#ifndef _ERROR_H_
#define _ERROR_H_

#include "Main.H"
#include <string>

class Errors
	{
	public:
		
		static void Warning(const std::string& a_string)
		{
			std::cout << "Warning: " << a_string << std::endl;
		}

		static void Abort(const std::string& a_string)
		{
			std::cout << "Fatal Error: " << a_string << std::endl;
			exit(1);
		}
		
		static void StepIn(const std::string& a_fname)
		{
			std::cout << a_fname << " -- step in" << std::endl;
		}
	};


#endif
