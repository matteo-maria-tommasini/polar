///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                 polar v1.0 - August the 15th 2019                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Copyright 2019 Matteo Tommasini                                          //
//                                                                           //
//  Licensed under the Apache License, Version 2.0 (the "License");          //
//  you may not use this file except in compliance with the License.         //
//  You may obtain a copy of the License at                                  //
//                                                                           //
//      http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                           //
//  Unless required by applicable law or agreed to in writing, software      //
//  distributed under the License is distributed on an "AS IS" BASIS,        //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
//  See the License for the specific language governing permissions and      //
//  limitations under the License.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef IOUTILS_H
#define IOUTILS_H

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

template<typename T>
std::vector<T> readitems(std::ifstream& infile, const int number_of_items)
{
   std::string line;
   std::vector<T> V;
   unsigned int to_be_read = number_of_items;

   std::getline(infile, line);
   while (to_be_read > 0)
   {
      boost::char_separator<char> sep(" ");
      boost::tokenizer< boost::char_separator<char> > tokens(line, sep);
      BOOST_FOREACH (const std::string& t, tokens) 
      { 
         T X;
         if (std::stringstream(t) >> X) { V.push_back(X); to_be_read -= 1; }      
      }
      if (to_be_read > 0) std::getline(infile, line);
      if (infile.eof()) 
      {
         std::cout << std::endl;
         std::cout << "Unexpected end of input file. Abort." << std::endl;
         exit (1);
      }
   }
   return(V);
}


std::string extract_ints(const std::ctype_base::mask category, 
                         const std::string& str, std::ctype<char> const& facet)
{
    using std::strlen;

    char const *begin = &str.front(),
               *end   = &str.back();

    auto res = facet.scan_is(category, begin, end);

    begin = &res[0];
    end   = &res[strlen(res)];

    return std::string(begin, end);
}

std::string extract_ints(const std::string& str)
{
    return extract_ints(std::ctype_base::digit, str,
         std::use_facet<std::ctype<char>>(std::locale("")));
}


#endif // IOUTILS_H

