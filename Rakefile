namespace :pkg do

  namespace :site do
    desc "Make the pkg website"
    task :build do
      sh 'Rscript -e "pkgdown::build_site(example = FALSE)"'
    end

    desc "Open the pkg website"
    task :open do
      sh 'xdg-open docs/index.html'
    end

    desc "Remove docs dir"
    task :clean do
      sh 'rm -rf docs/'
    end
  end

end
