# Intel 8051 Assembly Language Cheatsheet
---

## 🧮 Arithmetic Operations

### Addition
```asm
MOV A, #05H
MOV R0, #05H
ADD A, R0
```
- Adds 5 + 5, result in A.

### Subtraction
```asm
MOV A, #0AH
MOV R0, #05H
SUBB A, R0
```
- A = A - R0 - CY.

### Multiplication
```asm
MOV A, #05H
MOV B, #05H
MUL AB
```
- Result in A, overflow in B.

### Division
```asm
MOV A, #20H
MOV B, #05H
DIV AB
```
- Quotient in A, Remainder in B.

### Multi-byte Addition with Carry Check
```asm
MOV A, #92H
MOV R0, #23H
ADD A, R0
JNC L1
INC R7 ; Track carry
...
```
- `JNC`: Jump if No Carry.


---

## 🔗 Logical Operations / Digital Circuits
```asm
SETB ACC.0 ; A = 1
SETB ACC.1 ; B = 1
CLR ACC.2  ; C = 0
SETB ACC.3 ; D = 1

MOV C, ACC.0
ANL C, ACC.1
ORL C, ACC.2
CPL C
MOV ACC.7, C

MOV C, ACC.2
CPL C
ANL C, ACC.3
CPL C
ANL C, ACC.7 ; Final output F
```

---

## 🗃️ Stack Operations: Register Swapping
```asm
MOV R0, #25H
...
PUSH 0
PUSH 1
...
POP 0H
POP 1H
```
- Stack used to save and restore register values.
- SP increases on PUSH, decreases on POP.


---

## 🔁 ROM to RAM String Copy
```asm
MOV DPTR, #0200H
MOV R0, #40H
MOV R1, #0EH
LOOP: MOVC A,@A+DPTR
MOV @R0,A
INC DPTR
INC R0
DJNZ R1, LOOP
```
- `MOVC A,@A+DPTR`: Read from code memory.
- Copy string to internal RAM.

---

## 📌 Bit Addressable RAM
| Bit Addr | Byte Addr |
|----------|-----------|
| 42H      | 44H       |
| 67H      | 67H       |
| 0FH      | 01H       |
| 28H      | 28H       |
| 12H      | 12H       |
| 05H      | 00H       |

---

## 📤 Port Data Transfer with Delay
```asm
CPL P1.2
ACALL DELAY
...
DELAY: MOV R3,#200
HERE: DJNZ R3,HERE
RET
```
- Toggle P1.2 every ~436µs.

---

## ⏲️ Timer Programming (Square Wave)
```asm
MOV TMOD,#10H
MOV TL1,#06H
MOV TH1,#0FFH
SETB TR1
BACK: JNB TF1,BACK
CLR TR1
CPL P1.5
CLR TF1
SJMP AGAIN
```
- Timer-1, Mode-1: Generates 2kHz wave.

#### Question:
- How to generate square wave using Timer-1?

---

## ⏱️ Counter Programming (Pulse Count)
```asm
MOV TMOD,#01100000B ; Counter 1, Mode 2
MOV TH1,#0
SETB P3.5 ; T1 = input
SETB TR1
MOV A, TL1
MOV P2, A
```
- Counts pulses on pin T1 (P3.5), shows on Port 2.

---

## 📡 Serial Communication

### Transmit (Baud 9600)
```asm
MOV TMOD,#20H
MOV TH1,#-3 ;FDH
MOV SCON,#50H
SETB TR1
...
MOV SBUF, A
JNB TI, $ ; Wait
CLR TI
```

### Receive (Baud 4800)
```asm
MOV TMOD,#20H
MOV TH1,#-6
MOV SCON,#50H
SETB TR1
JNB RI, $ ; Wait
MOV A,SBUF
MOV P1,A
CLR RI
```

#### Question:
- Explain use of SCON and TH1 for baud rate.

---

## 🚨 Interrupt Programming (INT1)
```asm
ORG 0000H
LJMP MAIN
ORG 0013H ; ISR INT1
SETB P1.3 ; LED ON
DJNZ R3, $ ; delay
CLR P1.3 ; LED OFF
RETI
MAIN:
SETB TCON.3 ; INT1 edge-trigger
MOV IE,#10000100B ; Enable INT1
```
- INT1 (P3.3) triggers LED toggle via interrupt.
