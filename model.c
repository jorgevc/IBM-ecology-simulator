/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifdef SOI
	#define K(individual, Model) (Model->coagulation_factor)
#else
	#define K(individual, Model) (Model->coagulation_factor * (double)((individual.radio)*(individual.radio)))
#endif

//#define R(individual, Model) (individual.size_float/2.0)
//#define M(individual, Model) ((Model->Cm)/pow(2.0*Model->Cr,8.0/3.0)*pow(individual.size_float,2.6666))

#define R(individual, Model) (Model->Cr*pow(individual.size_float, 0.375))
#define M(individual, Model) (Model->Cm*individual.size_float)


