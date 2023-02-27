# demo.py

import ctypes
import llvmlite.binding as llvm

llvm.initialize()
llvm.initialize_native_target()
llvm.initialize_native_asmprinter()

strmod = """\
; ModuleID = "mod1"
target triple = "unknown-unknown-unknown"
target datalayout = ""

define double @"jit_func1"(double %"a")
{
entry:
  ;%".3" = fadd double 0x3ff0000000000000, %"a"
  %".3" = call double @llvm.cos.f64(double %"a")
  ret double %".3"
}

declare double    @llvm.cos.f64(double %Val)
"""

llmod = llvm.parse_assembly(strmod)

pmb = llvm.create_pass_manager_builder()
pmb.opt_level = 2
pass_manager = llvm.create_module_pass_manager()
pmb.populate(pass_manager)

pass_manager.run(llmod)

target_machine = llvm.Target.from_default_triple().create_target_machine()
exe_eng = llvm.create_mcjit_compiler(llmod, target_machine)
exe_eng.finalize_object()

fptr = exe_eng.get_function_address('jit_func1')

cfunc = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(fptr)

print(cfunc(4.0))
