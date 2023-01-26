
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
    else if fst(date) mod 4 =0 && fst(date) mod 400 =0 then
      Some 29
    else if fst(date) mod 4 =0 && fst(date) mod 100 !=0 then
      Some 29
    else
      Some 28

 (* let num = days_in_month(2100, 2);;  *)

let rec range (a:int) (b:int) (d1:int) (d2:int)  =
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
      
let wack = nth_day(2001, 165) 
;;


type country = {name:string; id:string; rates:float option list } 

(* let file = "csc330_a1.csv" *)

let read_file path =
let fp = open_in path in
let s = really_input_string fp (in_channel_length fp) in
close_in fp;
s



(* takes in a list of strings and converts each to floats *)
let to_float(nums: string list) = 
  let rec floater(numbers) =
    match numbers with
    | [] -> []
    | head :: tail -> Float.of_string_opt(head):: floater(tail)
  in floater(nums)


(* Takes in a list and outputs a country record *)
let rec sort(src:string list): country  = 
  match src with 
    | [] -> failwith("ERR: bad row")
    | _::[] -> failwith("ERR: bad row") 
    | head:: mid :: tail -> {name = head; id= mid; rates= to_float(tail)}


let get_records(file:string):country list =
  let output = read_file file in
  let new_list = String.split_on_char '\n' output in
  let rec loop(recs) =
    match recs with 
    | [] -> []
    | loopst::loopend -> sort(String.split_on_char ',' loopst) :: loop(loopend)
  in loop(new_list)


 let recs = get_records("csc330_a1.csv") 

;;

let check_rate(rate) =
  match rate with 
  | None -> 0
  | Some rate -> 1


let avail(data:country)=

let years = data.rates in
  let rec loop(numbers, sum) =
  match numbers with
  | [] -> sum
  | h :: t -> loop(t, sum +check_rate(h))
  in loop(years, 0)
  

  let af = {name = "Afghanistan"; id = "AFG";
  rates =
   [None; None; None; None; None; None; None; None; None; None; None; None;
    None; None; None; None; None; None; None; None; None; None; None; None;
    None; None; None; None; None; None; None; None; None; None; None; None;
    None; None; None; None; None; None; None; None; None; Some 12.68626872;
    Some 6.78459655; Some 8.680570785; Some 26.41866415; Some (-6.811161089);
    Some 2.178537524; Some 11.80418581; Some 6.441212809; Some 7.385771784;
    Some 4.673996035; Some (-0.661709165); Some 4.383891955;
    Some 4.975951506; Some 0.626149149; Some 2.302372515; None; None]}

let count = avail(af);;






let last(data:country)=
  let rec loop(rates, curr, year, curr_year) =
  match rates with
  | [] -> curr , year
  | h :: t -> (
    match h with
    | None -> loop(t, curr, year, curr_year+1)
    | Some float -> loop(t, h, curr_year+1, curr_year+1)
  )
  in loop(data.rates, Some 0.0, 1959, 1959)


let count = last(af);;




let minmax(record:country)=
  let rec loop(rates, min, max, year_min, year_max, curr_year)=
    match rates with 
    | [] -> (min,year_min),(max, year_max)
    | head :: tail ->
      (match head with 
      |None -> loop(tail, min, max, year_min, year_max, curr_year+1)
      |Some a ->(
        if a > max then 
          loop(tail, min, a, year_min, curr_year+1, curr_year+1) 
        else if a < min then
          loop(tail, a, max, curr_year+1, year_max, curr_year+1)
        else
          loop(tail, min, max, year_min, year_max, curr_year+1)
      )
      )in loop(record.rates, Float.max_float, Float.min_float, 1959,1959,1959)
    
let minny = minmax(af);;

let summarize(clist, identifier) =
  let rec loop(countries, iden) =
    match countries with 
    | [] -> None
    | hd::tl -> (
      if hd.id = iden then
        let latest = last(hd) in
        let str_last = Float.to_string(Option.get(fst(latest)))in
        let year_last = Int.to_string(snd(latest)) in
        let num_of =  Int.to_string(avail(hd)) in
        let min, max = minmax(hd) in
        let str_min = Float.to_string(fst(min))in
        let str_max = Float.to_string(fst(max))in
        let year_max = Int.to_string(snd(max)) in
        let year_min = Int.to_string(snd(min)) in
        let output = "Country: " ^ hd.name ^"(" ^hd.id^ ")\n" ^ 
        "Records available:" ^ num_of ^ "years\n" ^
        "Last record:" ^ year_last ^ "with rate of" ^ str_last ^"%\n" ^
        "Lowest rate:" ^ year_min ^ "with rate of" ^ str_min ^"%\n" ^
        "Highest rate:" ^ year_max ^ "with rate of" ^ str_max ^"%\n" 
      in Some output

      else 
        loop(tl, iden)  

    )  in loop(clist, identifier) 


let word = summarize(get_records("csc330_a1.csv"), "BOL");;

(* 
    let latest = last(hd) in
    let min, max = minmax(hd) in
    let str_min = Float.to_string(fst(min))in
    let str_max = Float.to_string(fst(max))in
    let year_max = Int.to_string(snd(max)) in
    let year_min = Int.to_string(snd(min)) in
    let output = "Country: " ^ hd.name ^'(' ^hd.id^ ")\n" ^  *)