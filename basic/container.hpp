#ifndef INCLUDEGUARD_CONTAINER_HPP
#define INCLUDEGUARD_CONTAINER_HPP

#include <ctime>     
#include <cstdlib>     
#include <iterator>
#include <functional>

#include "basic.hpp"



/** 
  * 
  * Container interface with on-the-fly implementation 
  * 
  * RETURNTYPE - data type of the values in the iteration 
  * ITERATOR_INTERNAL - data type of internal data for iterator
  * 
  */ 


template< class RETURNTYPE, class ITERATOR_INTERNAL>
class ContainerInterface
{
  
  private:
  
    /* internal members and function implementation - private */
    
    std::function<ITERATOR_INTERNAL()> m_initiate_iterator;
    
    std::function<void(ITERATOR_INTERNAL&)> m_destroy_iterator;
    
    std::function<bool(const ITERATOR_INTERNAL&)> m_is_iterator_finished;
    
    std::function<void(ITERATOR_INTERNAL&)> m_increment_iterator;
    
    std::function<RETURNTYPE(const ITERATOR_INTERNAL&)> m_dereference_iterator;
    
  public:
    
    /* constructor - public */
    /* simply initiates all the data */
    
    ContainerInterface(
        std::function<ITERATOR_INTERNAL()>                  a_initiate_iterator, 
        std::function<void(ITERATOR_INTERNAL&)>             a_destroy_iterator, 
        std::function<bool(const ITERATOR_INTERNAL&)>       a_is_iterator_finished, 
        std::function<void(ITERATOR_INTERNAL&)>             a_increment_iterator, 
        std::function<RETURNTYPE(const ITERATOR_INTERNAL&)> a_dereference_iterator
    ) : 
      m_initiate_iterator( a_initiate_iterator ),
      m_destroy_iterator ( a_destroy_iterator ),
      m_is_iterator_finished( a_is_iterator_finished ),
      m_increment_iterator( a_increment_iterator ),
      m_dereference_iterator( a_dereference_iterator )
    {}
    
    /* wrapper for function implementations - public */
    
    ITERATOR_INTERNAL initiate_iterator() const
    {
      return m_initiate_iterator();
    }
    
    void destroy_iterator(ITERATOR_INTERNAL& data ) const 
    {
      m_destroy_iterator( data );
    }
    
    bool is_iterator_finished(const ITERATOR_INTERNAL& data ) const 
    {
      return m_is_iterator_finished( data );
    }
    
    void increment_iterator(ITERATOR_INTERNAL& data ) const 
    {
      m_increment_iterator( data );
    }
    
    RETURNTYPE dereference_iterator(const ITERATOR_INTERNAL& data ) const 
    {
      return m_dereference_iterator( data );
    }
    
    
    /* member classes for iterators - public */
    
    struct range_end {};
    
    class iterator
    {
      
      private:
        
        const ContainerInterface& parentcontainer;
        
        ITERATOR_INTERNAL internaldata;
        
      public:
        
        iterator( const ContainerInterface& parentcontainer, ITERATOR_INTERNAL internaldata )
        : parentcontainer(parentcontainer), internaldata(internaldata)
        {};
        
        // dereference 
        const RETURNTYPE operator*() const
        {
          return parentcontainer.dereference_iterator( internaldata );
        }
        
        // increment prefix 
        const RETURNTYPE operator++()
        {
          parentcontainer.increment_iterator( internaldata );
          return parentcontainer.dereference_iterator( internaldata );
        }
        
        // increment postfix 
        const RETURNTYPE operator++(int)
        {
          RETURNTYPE&& tmp = parentcontainer.dereference_iterator( internaldata );
          parentcontainer.increment_iterator( internaldata );
          return tmp;
        }
        
        // check whether finished 
        bool operator!=( const range_end& ) const
        {
          return parentcontainer.is_iterator_finished( internaldata );
        }
        
    };
    
    /* for-each-loop interface - public */
    
    iterator begin() const { return iterator( *this, initiate_iterator() ); }
    
    const range_end end() const { return range_end(); }
    
};







#endif
