
.text
######## BARRETT ##########
.macro barrett a,kyb_Q,temp0,temp1 
    srli \temp0,\a,13         # u = a >> 13;//((uint32_t) a * sinv) >> 16;
    mul \temp1,\temp0,\kyb_Q   # u *= KYBER_Q;
    sub \a,\a,\temp1          # a -= u;   
.endm

.globl asm_barrett
.align 2
asm_barrett: 
    # Note the ABI: arguments a0-a7, return value is expected in a0 (and a1 if the case).      

    # args: a, Q
    barrett a0,a1,t0,t1    
    ret

