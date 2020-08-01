#ifndef INCLUDEGUARD_FIXEDARRAY_HPP
#define INCLUDEGUARD_FIXEDARRAY_HPP

#include <cassert>

#include <array>
#include <functional>
#include <utility>
#include <vector>


template<typename T>
class FixedArrayConstIterator;

template<typename T>
class FixedArrayIterator;




template <typename T>
class FixedArray 
{
    public:
    
        typedef T                   value_type;
        typedef unsigned int        size_type;
        typedef int                 difference_type;
        typedef value_type&         reference;
        typedef const value_type&   const_reference;
        
        typedef FixedArrayConstIterator<T>  const_iterator;
        typedef FixedArrayIterator<T>       iterator;
        
        static const int entrybytesize = sizeof(T); // ??? what if array

    private:
    
        size_type number_of_entries;
        T* data;
        
    public:
        
        explicit FixedArray( size_type length );
        explicit FixedArray( size_type length, T value );
        FixedArray( size_type length, std::function<T(size_type)> generator );
        FixedArray( size_type length, std::initializer_list<T> source );
//         FixedArray( size_type length, T* olddata );
  
        FixedArray( const FixedArray<T>& );
        FixedArray<T>& operator=( const FixedArray<T>& );
        FixedArray( FixedArray<T>&& );
        FixedArray<T>& operator=( FixedArray<T>&& );
        virtual ~FixedArray();
        
        const T& at( size_type i ) const;
        T& at( size_type i );
        
        const T& operator[]( size_type i ) const;
        T& operator[]( size_type i );
        
        size_type size() const;
        size_type length() const; // same as size()
        
        void swap( FixedArray<T>& other );

        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;
        const_iterator cbegin();
        const_iterator cend() const;
        
        void load( std::function<T(size_type)> );
        void load( const T* );
        
    public:
        
        
};





template<typename T>
FixedArray<T>::FixedArray( size_type length )
: number_of_entries(length),
  data(new T[length]())
{
    
}


template<typename T>
FixedArray<T>::FixedArray( size_type length, T value )
: number_of_entries(length),
  data( [length,value]()->T*{ T * foo = new T[length]; for( int i = 0; i < length; i++ ) foo[i] = value; return foo; }() )
{
    
}

template<typename T>
FixedArray<T>::FixedArray( size_type length, std::function<T(size_type)> generator )
: number_of_entries(length),
  data( [length,generator]()->T*{ T * foo = new T[length]; for( int i = 0; i < length; i++ ) foo[i] = generator(i); return foo; }() )
{
    
}

template<typename T>
FixedArray<T>::FixedArray( size_type length, std::initializer_list<T> source )
: number_of_entries(length),
  data(
    [length,source]()->T*{ 
      T * foo = new T[length]; 
      for( auto it = source.begin(); it != source.end(); it++ )
          foo[ it - source.begin() ] = *it; 
      return foo;
      }()
    )
{
    
}
// 
// template<typename T>
// FixedArray<T>::FixedArray( size_type length, T* olddata )
// : 
// {
//     
// }







template<typename T>
FixedArray<T>::FixedArray( const FixedArray<T>& other )
: FixedArray( other.length() )
{
    for( size_type i = 0; i < length(); i++ )
        data[i] = other.data[i];
}


template<typename T>
FixedArray<T>& FixedArray<T>::operator=( const FixedArray<T>& other )
{    
    if(this == &other)
      return *this;
    
    assert( length() == other.length() );
    
    if( length() == 0 ) 
        for( size_type i = 0; i < length(); i++ )
            data[i] = other.data[i];
    
    return *this;
}

template<typename T>
FixedArray<T>::FixedArray( FixedArray<T>&& other )
{
    number_of_entries = other.number_of_entries;
    data = other.data;
    
    other.number_of_entries = 0;
    delete[] other.data;
    other.data = nullptr;
}


template<typename T>
FixedArray<T>& FixedArray<T>::operator=( FixedArray<T>&& other )
{
    number_of_entries = 0;
    delete[] data;
    data = nullptr;
    
    number_of_entries = other.number_of_entries;
    data = other.data;
    
    other.number_of_entries = 0;
    delete[] other.data;
    other.data = nullptr;
    
    return *this;
}

template<typename T>
FixedArray<T>::~FixedArray(){
    number_of_entries = 0;
    delete[] data;
    data = nullptr;
}









template<typename T>
const T& FixedArray<T>::at( size_type i ) const
{
    assert( data != nullptr );
    assert( 0 <= i && i < length() );
    return data[i];
}

template<typename T>
T& FixedArray<T>::at( size_type i )
{
    assert( data != nullptr );
    assert( 0 <= i && i < length() );
    return data[i];
}
        
template<typename T>
const T& FixedArray<T>::operator[]( size_type i ) const
{
    assert( data != nullptr );
    assert( 0 <= i && i < length() );
    return data[i];
}

template<typename T>
T& FixedArray<T>::operator[]( size_type i )
{
    assert( data != nullptr );
    assert( 0 <= i && i < length() );
    return data[i];
}
        
template<typename T>
typename FixedArray<T>::size_type FixedArray<T>::size() const
{
    return length();
}

template<typename T>
typename FixedArray<T>::size_type FixedArray<T>::length() const
{
    return number_of_entries;
}
        
template<typename T>
void FixedArray<T>::swap( FixedArray<T>& other )
{
    auto temp_number_of_entries = other.number_of_entries;
    T* temp_data                = other.data;
    
    other.number_of_entries = number_of_entries;
    other.data              = data;
    
    number_of_entries = temp_number_of_entries;
    data              = temp_data;
}






template<typename T>
bool operator==( FixedArray<T>& one, FixedArray<T>& other )
{
    if( one.length() != other.length() )
        return false;
    
    for( typename FixedArray<T>::size_type i = 0; i < one.length(); i++ )
        if( one.at(i) != other.at(i) )
            return false;
    
    return true;
}

template<typename T>
bool operator!=( FixedArray<T>& one, FixedArray<T>& other )
{
    return !( one == other );
}


// Non member stuff:
// lexicographic comparision through < <= > >= 




template<typename T>
void swap( FixedArray<T>& one, FixedArray<T>& other )
{
    one.swap(other);
}






template<typename T>
void FixedArray<T>::load( std::function<T(size_type)> generator )
{
    for( typename FixedArray<T>::size_type i = 0; i < length(); i++ )
        data[i] = generator( i );
}

template<typename T>
void FixedArray<T>::load( const T* otherdata )
{
    for( typename FixedArray<T>::size_type i = 0; i < length(); i++ )
        data[i] = otherdata[i];
}




template<typename T>
FixedArrayIterator<T> FixedArray<T>::begin()
{
    return FixedArrayIterator<T>( this, 0 );
}

template<typename T>
FixedArrayIterator<T> FixedArray<T>::end() 
{
    return FixedArrayIterator<T>( this, length() );
}

template<typename T>
FixedArrayConstIterator<T> FixedArray<T>::begin() const
{
    return FixedArrayConstIterator<T>( this, 0 );
}

template<typename T>
FixedArrayConstIterator<T> FixedArray<T>::end() const
{
    return FixedArrayConstIterator<T>( this, length() );
}

template<typename T>
FixedArrayConstIterator<T> FixedArray<T>::cbegin()
{
    return begin();
}

template<typename T>
FixedArrayConstIterator<T> FixedArray<T>::cend() const
{
    return end();
}







template<typename Container, typename intT, typename extT>
class FixedArrayBaseIterator
{
    
    protected: 
        
        Container* parent;
        int index;
        
        FixedArrayBaseIterator( Container* parent, int index )
        : parent(parent), index(index) {}
        
    public:
        
        /* Comparison */
        
        bool operator!=( const FixedArrayBaseIterator<Container,intT,extT> other ) const
        {
            assert( this->parent == other.parent );
            return this->index != other.index;
        }
        
        bool operator==( const FixedArrayBaseIterator<Container,intT,extT> other ) const
        {
            return !( *this != other );
        }
        
        /* index operations */
        
        FixedArrayBaseIterator<Container,intT,extT>& operator++() 
        {
            index++;
            return *this;
        }
        
        FixedArrayBaseIterator<Container,intT,extT>& operator--() 
        {
            index--;
            return *this;
        }
        
        FixedArrayBaseIterator<Container,intT,extT>& operator++(int) 
        {
            index++;
            return FixedArrayBaseIterator<Container,intT,extT>( this->parent, index-1 );
        }
        
        FixedArrayBaseIterator<Container,intT,extT>& operator--(int) 
        {
            index--;
            return FixedArrayBaseIterator<Container,intT,extT>( this->parent, index+1 );
        }
        
//         FixedArrayBaseIterator<intT,extT>& operator+=( int i ) 
//         {
//             index+=i;
//             return *this;
//         }
//         
//         FixedArrayBaseIterator<intT,extT>& operator-=( int i ) 
//         {
//             index-=i;
//             return *this;
//         }
//         
//         FixedArrayBaseIterator<intT,extT>& operator+( int i ) 
//         {
//             auto ret = *this;
//             ret += i;
//             return i;
//         }
//         
//         FixedArrayBaseIterator<intT,extT>& operator-=( int i ) 
//         {
//             auto ret = *this;
//             ret += i;
//             return i;
//         }
        
        extT operator*() 
        {
            assert( 0 <= index && index < parent->length() );
            return parent->at(index);
        }
        
        extT& operator*() const
        {
            assert( 0 <= index && index < parent->length() );
            return parent->at(index);
        }
        
        
        
    
};

template<typename T>
class FixedArrayIterator
: public FixedArrayBaseIterator<FixedArray<T>,T,T>
{
    
    friend FixedArray<T>;
    
    protected:
    
        FixedArrayIterator( FixedArray<T>* parent, int index )
        : FixedArrayBaseIterator<FixedArray<T>,T,T>( parent, index )
        {}
    
};

template<typename T>
class FixedArrayConstIterator
: public FixedArrayBaseIterator<const FixedArray<T>,T,const T>
{
    
    friend FixedArray<T>;
    
    protected:
    
        FixedArrayConstIterator( const FixedArray<T>* parent, int index )
        : FixedArrayBaseIterator<const FixedArray<T>,T,const T>( parent, index )
        {}
    
};



// std::vector<const int> ggg(0);

// std::array<const int,3> goo = { 1,2 , 3};
// std::array<const int,0> hoo;



#endif
