
### An R-based web crawler / scraper 
### by Jack Davis
### jack.davis.sfu@gmail.com
### Last updated 2015-10-31
### For an explanation of the functions and parameters
### see: http://factotumjack.blogspot.ca/2015/10/an-r-based-web-crawler-scraper.html
 


library(XML)
library(stringr)
#library(RCurl)


######################################
########### Usage Examples

### Get all the within-site links from the front page of the National Post
links = get_clean_links("http://www.nationalpost.com/index.html", hasall="http://[a-zA-Z]+.(national|financial)post.com/")
head(links)

### Get the network structure found from exploring 5 web pages, starting at the front page of the National Post
url_seed= "http://www.nationalpost.com/index.html"
hasall = "http://[a-zA-Z]+.(national|financial)post.com/"
test = web_crawl(url_seed=url_seed, hasall=hasall, Npages=5, queue_method="Random")
names(test)
head(test$nodes)
head(test$edges)





#######################################
### Inner function to get the links from a single webpage, and prune out unwanted links

get_clean_links = function(input_url, prepend="", noindex=FALSE, nopost=TRUE, hasall=NULL,  hasnone=NULL, hassome=NULL)
{
require(XML)
require(stringr)
	
	links = as.character( try(getHTMLLinks(input_url), silent=TRUE)) ## Extract links
	### End early if we failed (can happen if a page doesn't exist)
	if(all(str_detect(links[1], "^Error"))){ print(links[1]); return("") }
	

	links = links[!str_detect(links,"pdf$")] # removes links to non-html documents
	links = links[!str_detect(links,"xml$")] # removes links to non-html documents

	if(noindex){	links = links[!str_detect(links,"index.html")]} # removes the links back to the main pages
	if(nopost){		links = links[!str_detect(links,"[ #?]")] } # removes links with POST information
	
	### Implement user defined filters
	if(length(hasall) > 0)
	{
		for(k in 1:length(hasall)){	links = links[ str_detect(links, hasall[k])]  }
	}
	if(length(hasnone) > 0)
	{
		for(k in 1:length(hasnone)){	links = links[ !str_detect(links, hasnone[k])]  }
	}
	if(length(hassome) > 0)
	{
		ingroup = numeric(0)
		for(k in 1:length(hassome)){	ingroup = union(ingroup, which( str_detect(links, hassome[k])))		}
		links = links[ingroup]
	}
	
	links = unique(links) # removes copies
	links = paste(prepend,links,sep="")
	
	return(links)
}



############################################
### Outer function - Crawls the web by scanning web pages for links and following those links.
### Returns two data frames, 'edges' and 'nodes', of the observed network of hyperlinks

web_crawl = function(url_seed="", Npages=1, queue_method = c("Random","BFS"), rand_seed=NULL, sleeptime=3, 
				hasall=NULL, hasnone=NULL, hassome=NULL, nopost=TRUE, noindex=FALSE, prepend="")
{
	if(length(rand_seed) == 1){set.seed(rand_seed)}
	count = 1
	failcount = 0
	
	if( !(queue_method[1] %in% c("BFS","Random"))){print("queue_method must be 'BFS' or 'Random'"); return(NULL)}
	if(queue_method[1] == "BFS"){	selection = 1 }
	if(queue_method[1] == "Random"){	selection = sample(1:length(url_seed),1)}
	this_url = url_seed[selection]
	link_queue = url_seed[-selection]
	
	pages_explored = c(this_url, rep("",Npages-1))
	degree_total = rep(0,Npages)
	degree_new = rep(0,Npages)
	all_edges_to = character(0)
	all_edges_from = character(0)

	### Loop until we get the requested number of pages or fail a bunch
	while(count < Npages & failcount < 10)
	{
		### Extract the links from the URL in focus
		new_links = get_clean_links(input_url=this_url, hasall=hasall, hasnone=hasnone, hassome=hassome, nopost=nopost, noindex=noindex, prepend="")
		new_links = new_links[new_links != ""] 
		degree_total[count] = length(new_links)
		
		### Update the vectors of edges
		all_edges_to = c(all_edges_to, new_links)
		all_edges_from = c(all_edges_from, rep(this_url, length(new_links)))
		
		### Remove previously explored links
		new_links = new_links[ !(new_links %in% pages_explored)]
		new_links = new_links[ !(new_links %in% link_queue)]
		degree_new[count] = length(new_links)
		
		
		### Update the vector of known links
		link_queue = c(link_queue, new_links)
		
		### Select a new focus URL
		if(queue_method[1] == "BFS"){	selection = 1 }
		if(queue_method[1] == "Random"){	selection = sample(1:length(link_queue),1)}
		this_url = link_queue[selection]
		link_queue = link_queue[-selection]
		
		
		link_queue = link_queue[!(link_queue %in% pages_explored)] ## Don't go in circles

		### Iteration recording / reporting
		count = count + 1
		pages_explored[count] = this_url
		print( paste(length(link_queue),   "queued,",count,":",this_url))
		
		### Rest for a time not to overload the server, collect garbage to clear RAM
		Sys.sleep(sleeptime)
		gc()
		if(length(link_queue) == 0){end} # end the loop prematurely if we run out of 
	}
	 
	
	### formatting node data and edge data
	nodes = data.frame(pages_explored, degree_total, degree_new)
	names(nodes)[1] = "address"
	nodes$samporder = 1:nrow(nodes)
	nodes$address = as.character(nodes$address)

	edges = data.frame(all_edges_from, all_edges_to)
	names(edges) = c("From","To")
	edges$From = as.character(edges$From)
	edges$To = as.character(edges$To)


	nodelist = c(nodes$address, setdiff(unique(edges$To),unique(edges$From)))

	### Assigning IDs to Nodes
	From_ID = rep(NA,nrow(edges))
	To_ID = rep(NA,nrow(edges))
	for(k in 1:length(nodelist))
	{
		target_url = nodelist[k]
		idx = which(edges$From == target_url)
		From_ID[idx] = k
		
		idx = which(edges$To == target_url)
		To_ID[idx] = k

	}

	edges$From_ID = From_ID
	edges$To_ID = To_ID

	return(list(nodes = nodes, edges = edges))
}

