
let fst3 (x, _, _) = x
let snd3 (_, x, _) = x
let thd3 (_, _, x) = x

let d1 = (2020, 8, 11)
let d2 = (2011, 11, 23)
let is_older((d1:int*int*int), (d2:int*int*int)): bool =
if fst3(d1) < fst3(d2) then 
    true
else if snd3(d1) < snd3(d2) && fst3(d1) = fst3(d2) then
    true
else if thd3(d1) < thd3(d2) && snd3(d1) = snd3(d2) && fst3(d1) = fst3(d2) then
    true
else 
    false


let oldest= is_older(d1,d2);;

let has31 = [1;3;5;7;8;10;12]
let has30 = [4;6;9;11] 
let days_in_month(date:int*int):int = 
    let month = snd(date) in
    if List.mem month has31 then
      31
    else if List.mem month has30 then
      30
    else if fst(date) mod 4 =0 || (fst(date) mod 400 =0 && fst(date) mod 100 =0) then
      29
    else
      28

(* let num = days_in_month(2100, 2);; *)

let rec range a b d1 d2  =
  if a > b then []
  else (d1, d2, a) :: range (a + 1 ) b d1 d2

let dates_in_month(date:int*int):(int*int*int) list =
let year  = fst(date) in
let month = snd(date) in  
let num = days_in_month(year,month) in
  

  let dates = range 1 num year month in dates



let blah = dates_in_month((2100,3));;







let rec sum((a:int), (date:int*int*int)):int  =
  if  a = snd3(date) then 0
  else
    days_in_month(fst3(date),a) + sum((a+1), date)


let num_of_days(date: int*int*int):int = 
  if snd3(date) = 1 then
    thd3(date)
  else
    sum(1, date) + thd3(date)


let bluh = num_of_days((2001, 1, 23));;
 



let rec neg((a:int), (n:int), (year:int)):int*int   =
  if  n <= 31 then a, n
  else a, (n - snd(neg((a+1), (n - days_in_month(year, (a+1))), year)) )


let nth_day((year:int), (n:int)):int*int*int = 
    if n <=31 then 
      (year, 1, n)
    else
      let month, day = neg(0, n, year)in 
      (year, month, day)
      



let wack = nth_day(2001, 50)
;;
   



  