// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__McIniData

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/nikikozyr/mcini/include/Nucleon.h"
#include "/home/nikikozyr/mcini/include/EventInitialState.h"
#include "/home/nikikozyr/mcini/include/URun.h"
#include "/home/nikikozyr/mcini/include/UEvent.h"
#include "/home/nikikozyr/mcini/include/UParticle.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Nucleon(void *p = 0);
   static void *newArray_Nucleon(Long_t size, void *p);
   static void delete_Nucleon(void *p);
   static void deleteArray_Nucleon(void *p);
   static void destruct_Nucleon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Nucleon*)
   {
      ::Nucleon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Nucleon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Nucleon", ::Nucleon::Class_Version(), "Nucleon.h", 23,
                  typeid(::Nucleon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Nucleon::Dictionary, isa_proxy, 4,
                  sizeof(::Nucleon) );
      instance.SetNew(&new_Nucleon);
      instance.SetNewArray(&newArray_Nucleon);
      instance.SetDelete(&delete_Nucleon);
      instance.SetDeleteArray(&deleteArray_Nucleon);
      instance.SetDestructor(&destruct_Nucleon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Nucleon*)
   {
      return GenerateInitInstanceLocal((::Nucleon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Nucleon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_EventInitialState(void *p = 0);
   static void *newArray_EventInitialState(Long_t size, void *p);
   static void delete_EventInitialState(void *p);
   static void deleteArray_EventInitialState(void *p);
   static void destruct_EventInitialState(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventInitialState*)
   {
      ::EventInitialState *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventInitialState >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EventInitialState", ::EventInitialState::Class_Version(), "EventInitialState.h", 16,
                  typeid(::EventInitialState), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventInitialState::Dictionary, isa_proxy, 4,
                  sizeof(::EventInitialState) );
      instance.SetNew(&new_EventInitialState);
      instance.SetNewArray(&newArray_EventInitialState);
      instance.SetDelete(&delete_EventInitialState);
      instance.SetDeleteArray(&deleteArray_EventInitialState);
      instance.SetDestructor(&destruct_EventInitialState);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventInitialState*)
   {
      return GenerateInitInstanceLocal((::EventInitialState*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EventInitialState*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_URun(void *p = 0);
   static void *newArray_URun(Long_t size, void *p);
   static void delete_URun(void *p);
   static void deleteArray_URun(void *p);
   static void destruct_URun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::URun*)
   {
      ::URun *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::URun >(0);
      static ::ROOT::TGenericClassInfo 
         instance("URun", ::URun::Class_Version(), "URun.h", 8,
                  typeid(::URun), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::URun::Dictionary, isa_proxy, 4,
                  sizeof(::URun) );
      instance.SetNew(&new_URun);
      instance.SetNewArray(&newArray_URun);
      instance.SetDelete(&delete_URun);
      instance.SetDeleteArray(&deleteArray_URun);
      instance.SetDestructor(&destruct_URun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::URun*)
   {
      return GenerateInitInstanceLocal((::URun*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::URun*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_UEvent(void *p = 0);
   static void *newArray_UEvent(Long_t size, void *p);
   static void delete_UEvent(void *p);
   static void deleteArray_UEvent(void *p);
   static void destruct_UEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::UEvent*)
   {
      ::UEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::UEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("UEvent", ::UEvent::Class_Version(), "UEvent.h", 12,
                  typeid(::UEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::UEvent::Dictionary, isa_proxy, 4,
                  sizeof(::UEvent) );
      instance.SetNew(&new_UEvent);
      instance.SetNewArray(&newArray_UEvent);
      instance.SetDelete(&delete_UEvent);
      instance.SetDeleteArray(&deleteArray_UEvent);
      instance.SetDestructor(&destruct_UEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::UEvent*)
   {
      return GenerateInitInstanceLocal((::UEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::UEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_UParticle(void *p = 0);
   static void *newArray_UParticle(Long_t size, void *p);
   static void delete_UParticle(void *p);
   static void deleteArray_UParticle(void *p);
   static void destruct_UParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::UParticle*)
   {
      ::UParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::UParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("UParticle", ::UParticle::Class_Version(), "UParticle.h", 10,
                  typeid(::UParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::UParticle::Dictionary, isa_proxy, 4,
                  sizeof(::UParticle) );
      instance.SetNew(&new_UParticle);
      instance.SetNewArray(&newArray_UParticle);
      instance.SetDelete(&delete_UParticle);
      instance.SetDeleteArray(&deleteArray_UParticle);
      instance.SetDestructor(&destruct_UParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::UParticle*)
   {
      return GenerateInitInstanceLocal((::UParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::UParticle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Nucleon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Nucleon::Class_Name()
{
   return "Nucleon";
}

//______________________________________________________________________________
const char *Nucleon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Nucleon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Nucleon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Nucleon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Nucleon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Nucleon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Nucleon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Nucleon*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr EventInitialState::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventInitialState::Class_Name()
{
   return "EventInitialState";
}

//______________________________________________________________________________
const char *EventInitialState::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventInitialState*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventInitialState::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventInitialState*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventInitialState::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventInitialState*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventInitialState::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventInitialState*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr URun::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *URun::Class_Name()
{
   return "URun";
}

//______________________________________________________________________________
const char *URun::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::URun*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int URun::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::URun*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *URun::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::URun*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *URun::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::URun*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr UEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *UEvent::Class_Name()
{
   return "UEvent";
}

//______________________________________________________________________________
const char *UEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int UEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *UEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *UEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr UParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *UParticle::Class_Name()
{
   return "UParticle";
}

//______________________________________________________________________________
const char *UParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int UParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *UParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *UParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Nucleon::Streamer(TBuffer &R__b)
{
   // Stream an object of class Nucleon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Nucleon::Class(),this);
   } else {
      R__b.WriteClassBuffer(Nucleon::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Nucleon(void *p) {
      return  p ? new(p) ::Nucleon : new ::Nucleon;
   }
   static void *newArray_Nucleon(Long_t nElements, void *p) {
      return p ? new(p) ::Nucleon[nElements] : new ::Nucleon[nElements];
   }
   // Wrapper around operator delete
   static void delete_Nucleon(void *p) {
      delete ((::Nucleon*)p);
   }
   static void deleteArray_Nucleon(void *p) {
      delete [] ((::Nucleon*)p);
   }
   static void destruct_Nucleon(void *p) {
      typedef ::Nucleon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Nucleon

//______________________________________________________________________________
void EventInitialState::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventInitialState.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EventInitialState::Class(),this);
   } else {
      R__b.WriteClassBuffer(EventInitialState::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventInitialState(void *p) {
      return  p ? new(p) ::EventInitialState : new ::EventInitialState;
   }
   static void *newArray_EventInitialState(Long_t nElements, void *p) {
      return p ? new(p) ::EventInitialState[nElements] : new ::EventInitialState[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventInitialState(void *p) {
      delete ((::EventInitialState*)p);
   }
   static void deleteArray_EventInitialState(void *p) {
      delete [] ((::EventInitialState*)p);
   }
   static void destruct_EventInitialState(void *p) {
      typedef ::EventInitialState current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EventInitialState

//______________________________________________________________________________
void URun::Streamer(TBuffer &R__b)
{
   // Stream an object of class URun.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(URun::Class(),this);
   } else {
      R__b.WriteClassBuffer(URun::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_URun(void *p) {
      return  p ? new(p) ::URun : new ::URun;
   }
   static void *newArray_URun(Long_t nElements, void *p) {
      return p ? new(p) ::URun[nElements] : new ::URun[nElements];
   }
   // Wrapper around operator delete
   static void delete_URun(void *p) {
      delete ((::URun*)p);
   }
   static void deleteArray_URun(void *p) {
      delete [] ((::URun*)p);
   }
   static void destruct_URun(void *p) {
      typedef ::URun current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::URun

//______________________________________________________________________________
void UEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class UEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(UEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(UEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_UEvent(void *p) {
      return  p ? new(p) ::UEvent : new ::UEvent;
   }
   static void *newArray_UEvent(Long_t nElements, void *p) {
      return p ? new(p) ::UEvent[nElements] : new ::UEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_UEvent(void *p) {
      delete ((::UEvent*)p);
   }
   static void deleteArray_UEvent(void *p) {
      delete [] ((::UEvent*)p);
   }
   static void destruct_UEvent(void *p) {
      typedef ::UEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::UEvent

//______________________________________________________________________________
void UParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class UParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(UParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(UParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_UParticle(void *p) {
      return  p ? new(p) ::UParticle : new ::UParticle;
   }
   static void *newArray_UParticle(Long_t nElements, void *p) {
      return p ? new(p) ::UParticle[nElements] : new ::UParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_UParticle(void *p) {
      delete ((::UParticle*)p);
   }
   static void deleteArray_UParticle(void *p) {
      delete [] ((::UParticle*)p);
   }
   static void destruct_UParticle(void *p) {
      typedef ::UParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::UParticle

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 216,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned short>*)0x0)->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete ((vector<unsigned short>*)p);
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] ((vector<unsigned short>*)p);
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

namespace ROOT {
   static TClass *vectorlENucleongR_Dictionary();
   static void vectorlENucleongR_TClassManip(TClass*);
   static void *new_vectorlENucleongR(void *p = 0);
   static void *newArray_vectorlENucleongR(Long_t size, void *p);
   static void delete_vectorlENucleongR(void *p);
   static void deleteArray_vectorlENucleongR(void *p);
   static void destruct_vectorlENucleongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Nucleon>*)
   {
      vector<Nucleon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Nucleon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Nucleon>", -2, "vector", 216,
                  typeid(vector<Nucleon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlENucleongR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Nucleon>) );
      instance.SetNew(&new_vectorlENucleongR);
      instance.SetNewArray(&newArray_vectorlENucleongR);
      instance.SetDelete(&delete_vectorlENucleongR);
      instance.SetDeleteArray(&deleteArray_vectorlENucleongR);
      instance.SetDestructor(&destruct_vectorlENucleongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Nucleon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Nucleon>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlENucleongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Nucleon>*)0x0)->GetClass();
      vectorlENucleongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlENucleongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlENucleongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Nucleon> : new vector<Nucleon>;
   }
   static void *newArray_vectorlENucleongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Nucleon>[nElements] : new vector<Nucleon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlENucleongR(void *p) {
      delete ((vector<Nucleon>*)p);
   }
   static void deleteArray_vectorlENucleongR(void *p) {
      delete [] ((vector<Nucleon>*)p);
   }
   static void destruct_vectorlENucleongR(void *p) {
      typedef vector<Nucleon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Nucleon>

namespace {
  void TriggerDictionaryInitialization_libMcIniData_Impl() {
    static const char* headers[] = {
"/home/nikikozyr/mcini/include/Nucleon.h",
"/home/nikikozyr/mcini/include/EventInitialState.h",
"/home/nikikozyr/mcini/include/URun.h",
"/home/nikikozyr/mcini/include/UEvent.h",
"/home/nikikozyr/mcini/include/UParticle.h",
0
    };
    static const char* includePaths[] = {
"/home/nikikozyr/Products/root/include",
"/home/nikikozyr/mcini/include",
"/home/nikikozyr/mcini",
"/home/nikikozyr/Products/root/include",
"/home/nikikozyr/mcini/build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libMcIniData dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$/home/nikikozyr/mcini/include/Nucleon.h")))  Nucleon;
class __attribute__((annotate("$clingAutoload$/home/nikikozyr/mcini/include/EventInitialState.h")))  EventInitialState;
class __attribute__((annotate("$clingAutoload$/home/nikikozyr/mcini/include/URun.h")))  URun;
class __attribute__((annotate("$clingAutoload$/home/nikikozyr/mcini/include/UEvent.h")))  UEvent;
class __attribute__((annotate("$clingAutoload$/home/nikikozyr/mcini/include/UParticle.h")))  UParticle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libMcIniData dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nikikozyr/mcini/include/Nucleon.h"
#include "/home/nikikozyr/mcini/include/EventInitialState.h"
#include "/home/nikikozyr/mcini/include/URun.h"
#include "/home/nikikozyr/mcini/include/UEvent.h"
#include "/home/nikikozyr/mcini/include/UParticle.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"EventInitialState", payloadCode, "@",
"Nucleon", payloadCode, "@",
"UEvent", payloadCode, "@",
"UParticle", payloadCode, "@",
"URun", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libMcIniData",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libMcIniData_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libMcIniData_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libMcIniData() {
  TriggerDictionaryInitialization_libMcIniData_Impl();
}
