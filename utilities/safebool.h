#ifndef __ITENSOR_SAFE_BOOL_H_
#define __ITENSOR_SAFE_BOOL_H_

namespace itensor {

#ifdef USE_CPP11

template <typename T>
class safe_bool
    {
    public:

    explicit operator bool() const
        {
        //std::cout << "Using C++11 safe bool idiom" << std::endl;
        return (static_cast<const T*>(this))->boolean_test();
        }

    protected:
    safe_bool() {}
    safe_bool(const safe_bool&) {}
    safe_bool& operator=(const safe_bool&) {return *this;}
    ~safe_bool() {}
    };


#else

template <typename T>
class safe_bool
    {
    protected:
    typedef void (safe_bool::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}

    public:
    operator bool_type() const 
        {
        //std::cout << "Using C++98 safe bool idiom" << std::endl;
        return (static_cast<const T*>(this))->boolean_test()
               ? &safe_bool::this_type_does_not_support_comparisons 
               : 0;
        }

    protected:
    safe_bool() {}
    safe_bool(const safe_bool&) {}
    safe_bool& operator=(const safe_bool&) {return *this;}
    ~safe_bool() {}
    };


template <typename T, typename U> 
bool
operator==(const safe_bool<T>& lhs, const safe_bool<U>& rhs) 
    {
    lhs.this_type_does_not_support_comparisons();   
    return false;
    }

template <typename T,typename U> 
bool
operator!=(const safe_bool<T>& lhs,const safe_bool<U>& rhs) 
    {
    lhs.this_type_does_not_support_comparisons();
    return false;   
    }

#endif

}; //namespace itensor
 

#endif //__ITENSOR_SAFE_BOOL_H_
