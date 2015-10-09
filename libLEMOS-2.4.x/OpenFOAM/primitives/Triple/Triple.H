/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Triple

Description
    An ordered triple of three objects of type \<T\> with first(), second()
    and third() elements.

\*---------------------------------------------------------------------------*/

#ifndef Triple_H
#define Triple_H

#include "FixedList.H"
#include "Istream.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Triple Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Triple
:
    public FixedList<Type, 3>
{

public:

    // Constructors

        //- Null constructor
        inline Triple()
        {}

        //- Construct from components
        inline Triple(const Type& f, const Type& s, const Type& t)
        {
            first() = f;
            second() = s;
            third() = t;
        }

        //- Construct from FixedList
        inline Triple(const FixedList<Type, 3>& fl)
        :
            FixedList<Type, 3>(fl)
        {}

        //- Construct from Istream
        inline Triple(Istream& is)
        :
            FixedList<Type, 3>(is)
        {}


    // Member Functions

        //- Return first
        inline const Type& first() const
        {
            return this->operator[](0);
        }

        //- Return first
        inline Type& first()
        {
            return this->operator[](0);
        }

        //- Return second
        inline const Type& second() const
        {
            return this->operator[](1);
        }

        //- Return second
        inline Type& second()
        {
            return this->operator[](1);
        }

        //- Return third
        inline const Type& third() const
        {
            return this->operator[](2);
        }

        //- Return third
        inline Type& third()
        {
            return this->operator[](2);
        }


        //- Return other
        inline const Pair<Type>& other(const Type& a) const
        {
            if (first() == second() || first() == third() || second() == third())
            {
                FatalErrorIn("Triple<Type>::other(const Type&) const")
                    << "Call to other only valid for Triple with differing"
                    << " elements:" << *this << abort(FatalError);
            }
            else if (first() == a)
            {
                return Pair<Type>(second(), third());
            }
            else if (second() == a)
            {
                return Pair<Type>(first(), third());
            }
            else
            {
                if (third() != a)
                {
                    FatalErrorIn("Triple<Type>::other(const Type&) const")
                        << "Triple " << *this
                        << " does not contain " << a << abort(FatalError);
                }
                return Pair<Type>(first(), second());
            }
        }


        //- Compare Triples
        //  Returning:
        //  -  0: different
        //  - +1: identical
        //  - -1: same pair, but reversed order
        static inline int compare(const Triple<Type>& a, const Triple<Type>& b)
        {
            if (a == b) 
            {
                return 1;
            }
            else if (a == reverse(b))
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
};


template<class Type>
Triple<Type> reverse(const Triple<Type>& p)
{
    return Triple<Type>(p.third(), p.second(), p.first());
}


template<class Type>
bool operator==(const Triple<Type>& a, const Triple<Type>& b)
{
    return (a.first() == b.first() && a.second() == b.second() && a.third() == b.third());
}


template<class Type>
bool operator!=(const Triple<Type>& a, const Triple<Type>& b)
{
    return !(a == b);
}


template<class Type>
bool operator<(const Triple<Type>& a, const Triple<Type>& b)
{
    return
    (
        (a.first() < b.first())
        && (a.second() < b.second())
     ||
        (
            !(b.first() < a.first())
         && !(b.second() < a.second())
         && a.third() < b.third()
        )
    );
}


template<class Type>
bool operator<=(const Triple<Type>& a, const Triple<Type>& b)
{
    return !(b < a);
}


template<class Type>
bool operator>(const Triple<Type>& a, const Triple<Type>& b)
{
    return (b < a);
}


template<class Type>
bool operator>=(const Triple<Type>& a, const Triple<Type>& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
