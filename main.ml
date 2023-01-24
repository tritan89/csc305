
let fst3 ((x:int), _, _) = x
let snd3 (_, x, _) = x
let thd3 (_, _, (x:int)) = x

let d1 = (2020, 8, 11)
let d2 = (2011, 11, 23)
let is_older((d1:int *int *int ), (d2:int *int *int )): bool =
if fst3(d1) < fst3(d2) then 
    true
else if snd3(d1) < snd3(d2) && fst3(d1) = fst3(d2) then
    true
else if thd3(d1) < thd3(d2) && snd3(d1) = snd3(d2) && fst3(d1) = fst3(d2) then
    true
else 
    false


(* let oldest= is_older(d1,d2);; *)

let has31 = [1;3;5;7;8;10;12]
let has30 = [4;6;9;11] 
let days_in_month(date:int*int ):int option   =
    if fst(date) > 3000 || fst(date) < 1 || snd(date) <1 || snd(date) >12  then 
      None
    else 
  
    let month = snd(date) in
    if List.mem month has31 then
      Some 31
    else if List.mem month has30 then
      Some 30
    else if fst(date) mod 4 =0 || (fst(date) mod 400 =0 && fst(date) mod 100 =0) then
      Some 29
    else
      Some 28

 (* let num = days_in_month(2100, 2);;  *)

let rec range (a:int  ) (b:int  ) (d1:int ) (d2:int )  =
  if a > b then []
  else (d1, d2, a) :: range (a + 1) b d1 d2

let dates_in_month(date:int *int ):(int *int *int ) list option  =
if fst(date) > 3000 || fst(date) < 1 || snd(date) <1 || snd(date) >12  then 
  None
else
let year  = fst(date) in
let month = snd(date) in  
let num = days_in_month(year,month) in
let clear = Option.get(num) in
let dates = range 1 clear year month in Some dates



(* let blah = dates_in_month((2100,3));; *)







let rec sum((a:int  ), (date:int *int *int )):int option  =

  if  a = snd3(date) then Some 0
  else
    let one  = Option.get( days_in_month(fst3(date),a)) in
    let two = Option.get((sum((a+1), date))) + one in
    Some two


let num_of_days(date: int *int *int ):int option   = 
if fst3(date) > 3000 || fst3(date) < 1 || snd3(date) <1 || snd3(date) >12 || thd3(date) <1 || thd3(date) >31 then 
  None
else
  if snd3(date) = 1 then
    Some (thd3(date))
  else
    let total = Option.get(sum(1, date)) in
    Some (total + thd3(date))


(* let bluh = num_of_days(2001, 11, 23);; *)
 



let rec neg((a:int ), (n:int ), (year:int )):int*int     =
  if  n <= 31 then  n, a
  else 
    let days = Option.get(days_in_month(year,a+1)) in
    neg((a+1), (n - days), year)


let nth_day((year:int ), (n:int )):int*int*int   = 
    if n <=31 then 
       (year, 1, n)
    else
      let day, month = neg(0, n, year) in 
       (year, month, day)
      
(* let wack = nth_day(2001, 165) *)
;;


type country = {id:string; name:string; rates:int list } ;;

let file = "csc330_a1.csv"

let read_file path =
let fp = open_in path in
let s = really_input_string fp (in_channel_length fp) in
close_in fp;
s




;;

let sort(src:string ):string

let get_records(file:string):country list =
  let output = read_file file in
  let clist = String.split_on_char ',' output in
    List.iter clist


