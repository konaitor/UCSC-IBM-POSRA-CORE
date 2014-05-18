#include "osra_polymer.h" 
#include <dirent.h>
using namespace std;
using namespace Magick;

string get_debug_path(string input_file) {
      time_t t;
      t = time(&t);
      struct tm *now = localtime(&t);
      stringstream ss;
      ss << (now->tm_mon + 1) << "-" << now->tm_mday << "_" << now->tm_hour << ":" << now->tm_min << "_";
      string file_name = string(basename((char*)input_file.c_str()));
      size_t extension = file_name.find_first_of(".");
      file_name.resize(extension);
      string debug_path = "./debug/" + ss.str() + file_name;
      if(!opendir(debug_path.c_str())) mkdir("./debug", 0777);
      if(!opendir(debug_path.c_str())) mkdir(debug_path.c_str(), 0777);
      string debug_name = debug_path + "/" + file_name;
      return debug_name;
}

void debug_log(string debug_log_name, double ave_bond_length, const vector<atom_t> atoms, const vector<bond_t> bonds){
      FILE *debug_file = fopen(debug_log_name.c_str(), "w");
      for(vector<bond_t>::const_iterator bond = bonds.begin(); bond != bonds.end(); ++bond) {
            if (bond->exists) {
                  fprintf(debug_file, "Bond %4d: ", (bond - bonds.begin())+1);
                  fprintf(debug_file, "[ %4d%3s%4d ] - [ %4d%3s%4d ]    ", bond->a+1, atoms[bond->a].label.c_str(), atoms[bond->a].anum,
                                                                           bond->b+1, atoms[bond->b].label.c_str(), atoms[bond->b].anum);
                  fprintf(debug_file, "[ x:%3d y:%3d ] - [ x:%3d y:%3d ]", (int)atoms[bond->a].x, (int)atoms[bond->a].y, (int)atoms[bond->b].x, (int)atoms[bond->b].y);
                  fprintf(debug_file, "\n");
            }
      }
      fprintf(debug_file, "\nAverage Bond Length: %f\n", ave_bond_length);
      fclose(debug_file);
}

void edit_smiles(string &s) {
      string tmp = s;
      vector<string> pieces;
      unsigned begin = 0;
      unsigned end = 0;
      while (begin < tmp.size() && end < tmp.size()) {
            end = tmp.find("[Po]", begin);
            if (end == string::npos) break;
            pieces.push_back(tmp.substr(begin, end - begin)); 
            begin = end + 4L;
      }
      begin = 0;
      end = 0;
      while (begin < tmp.size() && end < tmp.size()) {
            end = tmp.find("[Lv]", begin);
            if (end == string::npos) break;
            pieces.push_back(tmp.substr(begin, end - begin)); 
            begin = end + 4L;
      }
      int t = 0;
      for (vector<string>::iterator itor = pieces.begin(); itor != pieces.end(); ++itor)
            cout << string((t++%2==0)?"EG: ":"RU: ") << *itor << endl;
}

void  find_degree(Polymer &polymer, const vector<letters_t> letters, const vector<label_t> labels) {
      // Possible degree characters
      char possible_letters[] = "nmpxyNMPXY";
      // Somtimes OCRAD misinterprets a number for a character i.e. 50 -> 5O
      map<char, char> misread_numbers;
      misread_numbers['O'] = '0';
      misread_numbers['o'] = '0';
      // A collection of possible, or multiple degree characters, these could be a single character
      //  like n, or a string like "50"
      // The second in the pair, int, will represent the distance from the origin to associate to a
      //  a bracket
      vector<pair<string, double> > degrees;
      // First we check for strings, i.e. "50"
      for (vector<label_t>::const_iterator label = labels.begin(); label != labels.end(); ++label) {
            int degree;
            istringstream iss(label->a);
            iss >> degree;
            // Check if it's a valid number
            if (degree > 0) {
                  string s;
                  bool not_number = false;
                  for (string::const_iterator itor = label->a.begin(); itor != label->a.end(); ++itor) {
                        // Map characters to number characters
                        char c = misread_numbers[*itor];
                        // If c maps to something it will be nonzero, and we know it was a hit in our map
                        if (c) 
                              s.push_back(c); 
                        else {
                              int num;
                              istringstream iss(string(1L, *itor));
                              iss >> num;
                              if (num > 0)
                                    s.push_back(*itor);
                              else {
                                    not_number = true;
                                    break;
                              }
                        }
                  }
                  if (!not_number) {
                        double dis = sqrt((label->x1 * label->x1) + (label->y1 * label->y1));
                        degrees.push_back(make_pair(s, dis));
                  }
            }
      }
      // Next we check for characters, i.e. 'n'
      for (vector<letters_t>::const_iterator letter = letters.begin(); letter != letters.end(); ++letter) {
            if (letter->free) {
                  for (int i = 0; i < (sizeof(possible_letters) / sizeof(char)); ++i) {
                        // If a free (unassigned) letter is a possible degree specifier push it back!
                        if (letter->a == possible_letters[i]) {
                              double dis = sqrt((letter->x * letter->x) + (letter->y * letter->y));
                              degrees.push_back(make_pair(string(1L, letter->a), dis));
                        }
                  }
            }
      }
      // Debug print, need to formalize and associate with a bracket within a polymer
      for (vector<pair<string, double> >::iterator degree = degrees.begin(); degree != degrees.end(); ++degree) {
            cout << degree->first << endl;
      }
      // Right now it just takes the first degree that it could be associated with
      if (!degrees.empty()) {
            for (vector<pair<Bracket, Bracket> >::iterator bracket = polymer.brackets.begin(); bracket != polymer.brackets.end(); ++bracket) {
                  double brx = (double)bracket->second.get_bottom_right_x();
                  double bry = (double)bracket->second.get_bottom_right_y();
                  double b_dis = sqrt((brx * brx) + (bry * bry));
                  double threshold = (double)bracket->second.get_height();
                  for (vector<pair<string, double> >::iterator degree = degrees.begin(); degree != degrees.end(); ++degree) {
                        if (fabs(b_dis - degree->second) < threshold) {
                              bracket->first.set_degree(degree->first);
                              bracket->second.set_degree(degree->first);
                        }
                  }
            }
      }

}

void find_intersection(vector<bond_t> &bonds, const vector<atom_t> &atoms, vector<Bracket> &bracketboxes) {
      // As of now it only works on a single pair of brackets
      //if (bracketboxes.size() != 2) return;
      // Iterate through all of the bonds checking to see which ones intersect a detected Bracket
      for (vector<bond_t>::iterator bond = bonds.begin(); bond != bonds.end(); ++bond)
            if (bond->exists) {
                  for (int i = 0; i < bracketboxes.size(); i += 2) {
                        bool bracket0 = bracketboxes[i  ].intersects(*bond, atoms); // Check if the bonds intersects either bracket
                        bool bracket1 = bracketboxes[i+1].intersects(*bond, atoms);
                        //                  bool bracket0 = bracketboxes[0].intersects(*bond, atoms); // Check if the bonds intersects either bracket
                        //                  bool bracket1 = bracketboxes[1].intersects(*bond, atoms);
                        if (bracket0 || bracket1) {
                              bond->split = true;
                              // If the bond is split copy the corresponding orientation to the bond for later use
                              //if (bond->split) {
                              bond->bracket_orientation = (bracket0) ? bracketboxes[i].get_orientation() : bracketboxes[i+1].get_orientation();
                        }
                        //bond->bracket_orientation = (bracket0) ? bracketboxes[0].get_orientation() : bracketboxes[1].get_orientation();
                  }
            }
}

void pair_brackets(Polymer &polymer, const vector<Bracket> &brackets) {
      // i: Keeps track of bracket pairs, for every left there must be a right!
      int i = 0;
      for (vector<Bracket>::const_iterator bracket = brackets.begin(); bracket != brackets.end(); ++bracket) {
            if (i < 0) {
                  cerr << "Unmatched bracket" << endl;
                  break;
            }
            if (bracket->get_orientation() == 'l') {
                  polymer.set_polymer();
                  polymer.brackets.push_back(pair<Bracket, Bracket>(Bracket(), Bracket()));
                  polymer.brackets.back().first = *bracket;
                  ++i;
            } else {
                  --i;
                  polymer.brackets[i].second = *bracket;
            }
      }
      if (i > 0) {
            cerr << "Unmatched bracket" << endl;
      }
}

void  split_atom(vector<bond_t> &bonds, vector<atom_t> &atoms, int &n_atom, int &n_bond) {
      // Iterate through all of the bonds checking which ones have been marked to be split
      for(vector<bond_t>::iterator bond = bonds.begin(); bond != bonds.end(); ++bond) {
            if (bond->split) {
                  ++n_atom; // Increment the Total number of atoms, since were adding a polonium
                  ++n_bond; // Increment the total number of bonds, since were adding an atom
                  double x = (atoms[bond->a].x + atoms[bond->b].x) / 2; 
                  double y = (atoms[bond->a].y + atoms[bond->b].y) / 2;  
                  atom_t pseudo_atom(x, y, bond->curve); // Create a new atom
                  pseudo_atom.exists = true;
                  // Assign Polonium to left oriented brackets and Livermorium to right
                  pseudo_atom.label = (bond->bracket_orientation == 'l') ? "Po" : "Lv";
                  pseudo_atom.anum  = (bond->bracket_orientation == 'l') ? 84 : 116;
                  atoms.push_back(pseudo_atom);
                  // Create a new bond
                  bond_t newbond(atoms.size()-1, bond->b, bond->curve);
                  newbond.exists = true;
                  bonds.push_back(newbond);
                  bond->b = atoms.size() - 1;
            }
      }
}

void find_endpoints(Image detect, string debug_name, vector<pair<int, int> > &endpoints, int width, int height, vector<pair<pair<int, int>, pair<int, int> > > &bracketpoints) {
      const unsigned int SIDE_GROUP_SIZE = 2;
      const unsigned int BRACKET_MIN_SIZE = 5;
      for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                  int adj_groups = 0; // The number of side groups that contain at least 1 non-white pixel
                  vector<ColorGray> north, south, east, west;
                  vector<vector<ColorGray > > sides;

                  // Populate the side-groups vectors with the correct Color objects
                  for (int n = 1; n <= SIDE_GROUP_SIZE; n++ ) {
                        north.push_back (detect.pixelColor(i,j-n));
                        north.push_back (detect.pixelColor(i+1,j-n));
                        west.push_back  (detect.pixelColor(i-n,j));
                        west.push_back  (detect.pixelColor(i-n,j+1));
                        south.push_back (detect.pixelColor(i,j+(n+1)));
                        south.push_back (detect.pixelColor(i+1,j+(n+1)));
                        east.push_back  (detect.pixelColor(i+(n+1),j));
                        east.push_back  (detect.pixelColor(i+(n+1),j+1));
                  }
                  // Add each side-group vector the the sides vector
                  sides.push_back(north); sides.push_back(south); sides.push_back(east); sides.push_back(west);

                  // Check if each side group contains at least 1 non-white pixel
                  for (int c = 0; c < sides.size(); c++) {
                        for (int d = 0; d < sides.at(c).size(); d++) {
                              if (sides.at(c).at(d).shade() < 1) {
                                    adj_groups++;
                                    break;
                              }
                        }
                  }
                  if (adj_groups == 0) continue;

                  // If the current pixel, or the pixel to its immediate right are non-white,
                  // and 3 of the adjacent groups are completely white, then add the current
                  // pixel (or immediate right) to list of endpoints.
                  ColorGray current_pixel = detect.pixelColor(i,j);
                  ColorGray current_right = detect.pixelColor(i+1,j);
                  if (adj_groups == 1 && (current_pixel.shade() < 1 || current_right.shade() < 1)) {
                        if (current_pixel.shade() < 1 && current_right.shade() == 1)
                              endpoints.push_back(make_pair(i,j));
                        else if (current_pixel.shade() == 1 && current_right.shade() < 1)
                              endpoints.push_back(make_pair(i+1,j));
                  }
            }
      }

      // For each endpoint, loop through all other endpoints
      for (int i = 0; i < endpoints.size(); i++) {
            // Uncomment next line to show all endpoints:
            detect.pixelColor (endpoints.at(i).first, endpoints.at(i).second, "red");
            for (int j = 0; j < endpoints.size(); j++) {
                  // If x values are equivalent (+/- 1) and the endpoints are reasonable distance apart
                  // then procede to next test
                  if ((endpoints.at(i).first == endpoints.at(j).first ||
                                    endpoints.at(i).first == endpoints.at(j).first + 1 ||
                                    endpoints.at(i).first == endpoints.at(j).first - 1) && i != j &&
                              abs (endpoints.at(i).second - endpoints.at(j).second) > BRACKET_MIN_SIZE) {
                        const unsigned int lower_index = (endpoints.at(i).second < endpoints.at(j).second ? i : j);
                        const unsigned int upper_index = (endpoints.at(i).second > endpoints.at(j).second ? i : j);
                        const unsigned int mid_y = (endpoints.at(i).second + endpoints.at(j).second) / 2;
                        const unsigned int qtr0_y = (mid_y + endpoints.at(lower_index).second)/2;
                        const unsigned int qtr1_y = (mid_y + endpoints.at(upper_index).second)/2;
                        ColorGray qtr0 = detect.pixelColor(endpoints.at(upper_index).first, qtr0_y);
                        ColorGray qtr1 = detect.pixelColor(endpoints.at(lower_index).first, qtr0_y);
                        ColorGray qtr2 = detect.pixelColor(endpoints.at(upper_index).first, qtr1_y);
                        ColorGray qtr3 = detect.pixelColor(endpoints.at(lower_index).first, qtr1_y);

                        unsigned int non_white = 0;
                        bool intersects = false;
                        for (unsigned int k = qtr0_y; k <= qtr1_y; k++) {
                            ColorGray c = detect.pixelColor(endpoints.at(i).first, k);
                            if (c.shade() < 1) {
                                intersects = true;
                                break;
                            }
                        }
                        for (unsigned int k = endpoints.at(lower_index).second; k < endpoints.at(upper_index).second; k++) {
                            ColorGray c = detect.pixelColor(endpoints.at(i).first, k);
                            if (c.shade() < 1) non_white++;
                        }
                        // If at least one of the middle pixels are non-white and both of the quarter
                        // pixels are white, then add pair as endpoint
                        if (intersects && non_white < 10 && qtr0.shade() == 1 && qtr1.shade() == 1 && qtr2.shade() == 1 && qtr3.shade() == 1) {
                              // Drawing the green line between endpoints
                              for (int b = endpoints.at(lower_index).second; b < endpoints.at(upper_index).second; b++)
                                    detect.pixelColor (endpoints.at(j).first, b, "green");
                              detect.pixelColor (endpoints.at(i).first, endpoints.at(i).second, "blue");
                              detect.pixelColor (endpoints.at(j).first, endpoints.at(j).second, "blue");
                              // Printing coordinates of bracket
                              //cout << endpoints.at(i).first << ", " << endpoints.at(i).second << endl;
                              //cout << endpoints.at(j).first << ", " << endpoints.at(j).second << endl;
                              bracketpoints.push_back(make_pair(endpoints.at(i), endpoints.at(j)));
                                    
                              endpoints.erase(endpoints.begin() + (i > j ? i : j));
                              endpoints.erase(endpoints.begin() + (i < j ? i : j));
                        }                                   
                  }
            }
      }
      detect.write(debug_name + "_endpoints_detect.gif");
}

void find_brackets(Image &img, string debug_name, vector<Bracket> &bracketboxes) { 
      vector<pair<int, int> > endpoints;
      vector<pair<pair<int, int>,pair<int, int> > > bracketpoints;
      // Find endpoints in the image
      find_endpoints(img, debug_name, endpoints, img.columns(), img.rows(), bracketpoints);
      //if(bracketpoints.size() != 2) return;
      // Iterate over endpoints and convert them into Brackets
      for(vector<pair<pair<int, int>, pair<int, int> > >::iterator itor = bracketpoints.begin(); itor != bracketpoints.end(); ++itor)
            bracketboxes.push_back(Bracket(itor->first, itor->second, img)); 

      // Remove the brackets from the image entirely
      //bracketboxes[0].remove_brackets(img);
      //bracketboxes[1].remove_brackets(img);
      for (int i = 0; i < bracketboxes.size(); i++) {
            bracketboxes[i].remove_brackets(img);
      }

}

void plot_points(Image &img, const vector<point> &points) {
      for(vector<point>::const_iterator itor = points.begin(); itor != points.end(); ++itor)
            if(itor->color != CLEAR) {
                  int x = itor->x, y = itor->y;
                  if(x < 0)             x = 0;
                  if(x > img.columns()) x = img.columns() - 1;
                  if(y < 0)             y = 0;
                  if(y > img.rows())    y = img.rows() - 1;
                  img.pixelColor(x, y, colors[itor->color]);
            }
}

void plot_atoms(Image &img, const vector<atom_t> &atoms, const std::string color) {
      for(vector<atom_t>::const_iterator itor = atoms.begin(); itor != atoms.end(); ++itor)
            if(itor->exists) img.pixelColor(itor->x, itor->y, color);
}

void plot_bonds(Image &img, const vector<bond_t> &bonds, const vector<atom_t> atoms, const std::string color) {
      for(vector<bond_t>::const_iterator itor = bonds.begin(); itor != bonds.end(); ++itor)
            if(itor->exists) {
                  int x = (atoms[itor->a].x + atoms[itor->b].x) / 2;
                  int y = (atoms[itor->a].y + atoms[itor->b].y) / 2;
                  if(itor->split) img.pixelColor(x, y, "green");
                  else img.pixelColor(x, y, color);
            }
}

void plot_letters(Image &img, const vector<letters_t> &letters, const std::string color) {
      for(vector<letters_t>::const_iterator itor = letters.begin(); itor != letters.end(); ++itor)
            img.pixelColor(itor->x, itor->y, color);
}

void plot_labels(Image &img, const vector<label_t> &labels, const std::string color) {
      for(vector<label_t>::const_iterator itor = labels.begin(); itor != labels.end(); ++itor) {
            img.pixelColor(itor->x1, itor->y1, color);
            img.pixelColor(itor->x2, itor->y2, color);
      }
}

void plot_all(Image img, string debug_name, const int boxn, const string id, const vector<atom_t> atoms, const vector<bond_t> bonds, const vector<letters_t> letters, const vector<label_t> labels) {
      plot_atoms(img, atoms, "blue");
      plot_bonds(img, bonds, atoms, "red");
      plot_letters(img, letters, "orange");
      plot_labels(img, labels, "pink");
      ostringstream ss;
      ss << boxn;
      img.write(debug_name + "_atoms_bonds_labels_" + ss.str() + ".gif");
}
