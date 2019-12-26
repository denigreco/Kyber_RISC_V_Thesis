# extern void ntt_merged(uint16_t *p,uint16_t *zetas, uint32_t KQ, uint32_t KQ4,uint32_t qinv,new_zeta);

.text

######## MONTGOMERRY #########
#return value in arg1
.macro montgomery a,qinv,kyb_Q,temp0,temp1,mont_const
  mul \temp0,\a,\qinv        #u = (a * qinv);   
                              #first 32 bits of the multiplication in temp0  
  and \temp0,\temp0,\mont_const    #u &= ((1<<rlog)-1);    
  mul \temp1,\temp0,\kyb_Q    #u *= KYBER_Q;  
  add \a,\a,\temp1            #a = a + u;
  srli \a,\a,18               #return a >> rlog; (uint16_t)  
.endm

######## BARRETT ##########
.macro barrett a,kyb_Q,temp0,temp1 
    srli \temp0,\a,13         # u = a >> 13;//((uint32_t) a * sinv) >> 16;
    mul \temp1,\temp0,\kyb_Q   # u *= KYBER_Q;
    sub \a,\a,\temp1          # a -= u;   
.endm

####### BUTTERFLY ODD ##########
.macro butterfly_odd left,right,zeta,q4,qinv,kyb_Q,temp0,temp1,t,mont_const   
    #left  = p[i + aux]
    #right  = p[i + aux + mod]

    mul \t,\zeta,\right       
    montgomery \t,\qinv,\kyb_Q,\temp0,\temp1,\mont_const   #t = montgomery_reduce((uint32_t)zeta * p[i + aux + mod]);   
                         
    #right term before reduction: p[i + aux] + 4*KYBER_Q - t    
    add \right,\left,\q4     #p[i + aux] + 4*KYBER_Q
    sub \right,\right,\t      #p[i + aux] + 4*KYBER_Q - t                        
           
    #left term before reduction
    add \left,\left,\t

    #barret reduction right
    barrett \right,\kyb_Q,\temp0,\temp1
                
    # barret reduction left - even level only!
    #barrett \left,\kyb_Q,\temp0,\temp1      
     
.endm


####### BUTTERFLY EVEN ##########
.macro butterfly_even left,right,zeta,q4,qinv,kyb_Q,temp0,temp1,t,mont_const   
    #left  = p[i + aux]
    #right  = p[i + aux + mod]

    mul \t,\zeta,\right       
    montgomery \t,\qinv,\kyb_Q,\temp0,\temp1,\mont_const   #t = montgomery_reduce((uint32_t)zeta * p[i + aux + mod]);   
                         
    #right term before reduction: p[i + aux] + 4*KYBER_Q - t    
    add \right,\left,\q4     #p[i + aux] + 4*KYBER_Q
    sub \right,\right,\t      #p[i + aux] + 4*KYBER_Q - t                        
           
    #left term before reduction
    add \left,\left,\t

    #barret reduction right
    barrett \right,\kyb_Q,\temp0,\temp1
                
    # barret reduction left - even level only!
    barrett \left,\kyb_Q,\temp0,\temp1      
     
.endm




.globl ntt_merged
.align 2
ntt_merged: #args (p,zetas, KYBER_Q, KYBER_Q<<2, 7679U)
    # Note the ABI: arguments a0-a7, return value is expected in a0 (and a1 if the case). 
    # If used s0-s11, ra,sp - save to stack and restore after 
    # If the values of t0-t6 needed after calling a fnct, save it before    

    addi sp, sp, -52 	# Make space for 13 words on the stack
    sw s0,0(sp)
    sw s1,4(sp)
    sw s2,8(sp)
    sw s3,12(sp)
    sw s4,16(sp)
    sw s5,20(sp)
    sw s6,24(sp)
    sw s7,28(sp)
    sw s8,32(sp)
    sw s9,36(sp)
    sw s10,40(sp)
    sw s11,44(sp)
    sw ra,48(sp)


 
    # ra - montgomery constant 
    addi ra,x0,1          
    slli ra,ra,18
    addi ra,ra,-1       # ra = (1<<18)-1 


##############  LEVEL 7_6_5_4  ##############
    addi a0,a0,32   #p=p[16]
    addi t6,x0,15

do:     
    addi a0,a0,-2   #p-- 


    #load 16 butterflies values 
    lhu s0,0(a0)        # s0=p[0]
    lhu s1,32(a0)        # s1=p[16]
    lhu s2,64(a0)        # s2=p[32]
    lhu s3,96(a0)        # s3=p[48]
    lhu s4,128(a0)        # s4=p[64]
    lhu s5,160(a0)        # s5=p[80]
    lhu s6,192(a0)        # s6=p[96]
    lhu s7,224(a0)        # s7=p[112]
    lhu s8,256(a0)        # s8=p[128]
    lhu s9,288(a0)        # s9=p[144]
    lhu s10,320(a0)        # s10=p[160]
    lhu s11,352(a0)        # s11=p[176]
    lhu t5,384(a0)        # t5=p[192]
    lhu a6,416(a0)        # a6=p[208]
    lhu a7,448(a0)        # a7=p[224]
    lhu t3,480(a0)        # t3=p[240]


    #########  LEVEL 7  ###########
    lhu t4,2(a1)     #load zeta
    butterfly_odd s0,s8,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s1,s9,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s2,s10,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s3,s11,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s4,t5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s5,a6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s6,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s7,t3,t4,a3,a4,a2,t0,t1,t2,ra
    #########  LEVEL 6  ###########
    lhu t4,4(a1)     #load zeta
    butterfly_even s0,s4,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s1,s5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s2,s6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s3,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,6(a1)     #load zeta
    butterfly_even s8,t5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s9,a6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s10,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s11,t3,t4,a3,a4,a2,t0,t1,t2,ra
    #########  LEVEL 5  ###########
    lhu t4,8(a1)     #load zeta
    butterfly_odd s0,s2,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s1,s3,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,10(a1)     #load zeta
    butterfly_odd s4,s6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s5,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,12(a1)     #load zeta
    butterfly_odd s8,s10,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s9,s11,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,14(a1)     #load zeta
    butterfly_odd t5,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd a6,t3,t4,a3,a4,a2,t0,t1,t2,ra
    #########  LEVEL 4  ###########
    lhu t4,16(a1)     #load zeta
    butterfly_even s0,s1,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,18(a1)     #load zeta
    butterfly_even s2,s3,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,20(a1)     #load zeta
    butterfly_even s4,s5,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,22(a1)     #load zeta
    butterfly_even s6,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,24(a1)     #load zeta
    butterfly_even s8,s9,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,26(a1)     #load zeta
    butterfly_even s10,s11,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,28(a1)     #load zeta
    butterfly_even t5,a6,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,30(a1)     #load zeta
    butterfly_even a7,t3,t4,a3,a4,a2,t0,t1,t2,ra

    #store 16 butterflies values 
    sh s0,0(a0)        # p[0]=s0
    sh s1,32(a0)        # p[16]=s1
    sh s2,64(a0)        # p[32]=s2
    sh s3,96(a0)        # p[48]=s3
    sh s4,128(a0)        # p[64]=s4
    sh s5,160(a0)        # p[80]=s5
    sh s6,192(a0)        # p[96]=s6
    sh s7,224(a0)        # p[112]=s7
    sh s8,256(a0)        # p[128]=s8
    sh s9,288(a0)        # p[144]=s9
    sh s10,320(a0)        # p[160]=s10
    sh s11,352(a0)        # p[176]=s11
    sh t5,384(a0)        # p[192]=t5
    sh a6,416(a0)        # p[208]=a6
    sh a7,448(a0)        # p[224]=a7
    sh t3,480(a0)        # p[240]=t3


    addi t6,t6,-1
while:
    bge t6,zero,do     #loop while t6>=0 //16 loops


##############  END LEVEL 7_6_5_4  ##############




##############  LEVEL 3-0  ##############
    addi t6,x0,16

loop_16:
    #load 16 butterflies values 
    lhu s0,0(a0)        # s0=p[0]
    lhu s1,2(a0)        # s1=p[1]
    lhu s2,4(a0)        # s2=p[2]
    lhu s3,6(a0)        # s3=p[3]
    lhu s4,8(a0)        # s4=p[4]
    lhu s5,10(a0)        # s5=p[5]
    lhu s6,12(a0)        # s6=p[6]
    lhu s7,14(a0)        # s7=p[7]
    lhu s8,16(a0)        # s8=p[8]
    lhu s9,18(a0)        # s9=p[9]
    lhu s10,20(a0)        # s10=p[10]
    lhu s11,22(a0)        # s11=p[11]
    lhu t5,24(a0)        # t5=p[12]
    lhu a6,26(a0)        # a6=p[13]
    lhu a7,28(a0)        # a7=p[14]
    lhu t3,30(a0)        # t3=p[15]


##############  LEVEL 3  ##############
    lhu t4,0(a5)     #load zeta new
    butterfly_odd s0,s8,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s1,s9,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s2,s10,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s3,s11,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s4,t5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s5,a6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s6,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s7,t3,t4,a3,a4,a2,t0,t1,t2,ra
##############  LEVEL 2  ##############
    lhu t4,2(a5)     #load zeta new
    butterfly_even s0,s4,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s1,s5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s2,s6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s3,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,4(a5)     #load zeta new
    butterfly_even s8,t5,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s9,a6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s10,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_even s11,t3,t4,a3,a4,a2,t0,t1,t2,ra
##############  LEVEL 1  ##############
    lhu t4,6(a5)     #load zeta new
    butterfly_odd s0,s2,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s1,s3,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,8(a5)     #load zeta new
    butterfly_odd s4,s6,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s5,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,10(a5)     #load zeta new
    butterfly_odd s8,s10,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd s9,s11,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,12(a5)     #load zeta new
    butterfly_odd t5,a7,t4,a3,a4,a2,t0,t1,t2,ra
    butterfly_odd a6,t3,t4,a3,a4,a2,t0,t1,t2,ra
##############  LEVEL 0  ##############
    lhu t4,14(a5)     #load zeta new
    butterfly_even s0,s1,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,16(a5)     #load zeta new
    butterfly_even s2,s3,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,18(a5)     #load zeta new
    butterfly_even s4,s5,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,20(a5)     #load zeta new
    butterfly_even s6,s7,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,22(a5)     #load zeta new
    butterfly_even s8,s9,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,24(a5)     #load zeta new
    butterfly_even s10,s11,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,26(a5)     #load zeta new
    butterfly_even t5,a6,t4,a3,a4,a2,t0,t1,t2,ra
    lhu t4,28(a5)     #load zeta new
    butterfly_even a7,t3,t4,a3,a4,a2,t0,t1,t2,ra

    #store 16 butterflies values 
    sh s0,0(a0)        # p[0]=s0
    sh s1,2(a0)        # p[1]=s1
    sh s2,4(a0)        # p[2]=s2
    sh s3,6(a0)        # p[3]=s3
    sh s4,8(a0)        # p[4]=s4
    sh s5,10(a0)        # p[5]=s5
    sh s6,12(a0)        # p[6]=s6
    sh s7,14(a0)        # p[7]=s7
    sh s8,16(a0)        # p[8]=s8
    sh s9,18(a0)        # p[9]=s9
    sh s10,20(a0)        # p[10]=s10
    sh s11,22(a0)        # p[11]=s11
    sh t5,24(a0)        # p[12]=t5
    sh a6,26(a0)        # p[13]=a6
    sh a7,28(a0)        # p[14]=a7
    sh t3,30(a0)        # p[15]=t3


    addi a0,a0,32      #p+=16
    addi a5,a5,30
    addi t6,t6,-1
    bne t6,x0,loop_16  #jump back if >=0
##############  END LEVEL 3-0  ##############
    lw s0,0(sp)
    lw s1,4(sp)
    lw s2,8(sp)
    lw s3,12(sp)
    lw s4,16(sp)
    lw s5,20(sp)
    lw s6,24(sp)
    lw s7,28(sp)
    lw s8,32(sp)
    lw s9,36(sp)
    lw s10,40(sp)
    lw s11,44(sp)
    lw ra,48(sp)
    addi sp, sp, 52 	# Clear 13 words from the stack


    ret
